import os
import binascii
import unittest


from x25519.diffie_hellman import EllipticCurveDiffieHellman
from nacl.bindings import crypto_scalarmult, crypto_box_beforenm
from nacl.public import PrivateKey
from Crypto.Cipher import Salsa20

class TestDiffieHellmanVectors(unittest.TestCase):
    def setUp(self):
        self.base_point = b'\x09' + (b'\x00' * 31)

    def test_alice_public_key(self):
        # Test vector from RFC 7748 / RFC 8439:
        # Alice's private key, a:
        alice_private_hex = "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a"
        # Expected Alice's public key, X25519(a, 9):
        expected_alice_public_hex = "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a"
        
        alice_private = binascii.unhexlify(alice_private_hex)
        expected_alice_public = binascii.unhexlify(expected_alice_public_hex)
        
        alice = EllipticCurveDiffieHellman(private_key=alice_private)
        self.assertEqual(alice.public_key, expected_alice_public,
                         msg="Alice public key does not match the test vector.")

    def test_bob_public_key(self):
        # Bob's private key, b:
        bob_private_hex = "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb"
        # Expected Bob's public key, X25519(b, 9):
        expected_bob_public_hex = "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f"
        
        bob_private = binascii.unhexlify(bob_private_hex)
        expected_bob_public = binascii.unhexlify(expected_bob_public_hex)
        
        bob = EllipticCurveDiffieHellman(private_key=bob_private)
        self.assertEqual(bob.public_key, expected_bob_public,
                         msg="Bob public key does not match the test vector.")

    def test_shared_secret(self):
        # Using the same test vectors as above:
        alice_private_hex = "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a"
        bob_private_hex = "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb"
        # Expected shared secret, K:
        expected_shared_hex = "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742"
        
        alice_private = binascii.unhexlify(alice_private_hex)
        bob_private = binascii.unhexlify(bob_private_hex)
        expected_shared = binascii.unhexlify(expected_shared_hex)
        
        alice = EllipticCurveDiffieHellman(private_key=alice_private)
        bob = EllipticCurveDiffieHellman(private_key=bob_private)
        
        # Alice computes the shared secret using Bob's public key.
        alice_shared = alice.generate_shared_secret(bob.public_key)
        # Bob computes the shared secret using Alice's public key.
        bob_shared = bob.generate_shared_secret(alice.public_key)
        
        self.assertEqual(alice_shared, expected_shared,
                         msg="Alice's computed shared secret does not match the test vector.")
        self.assertEqual(bob_shared, expected_shared,
                         msg="Bob's computed shared secret does not match the test vector.")
        self.assertEqual(alice_shared, bob_shared,
                         msg="Alice and Bob's shared secrets do not match.")

    def test_compare_with_pynacl_public_key(self):
        # Generate a random private key.
        private = os.urandom(32)
        
        # Compute public key using our DiffieHellman abstraction.
        my_dh = EllipticCurveDiffieHellman(private_key=private)
        my_public = my_dh.public_key
        
        # Compute public key using PyNaCl's crypto_scalarmult.
        py_public = crypto_scalarmult(private, self.base_point)
        
        self.assertEqual(my_public, py_public,
                         msg="Public key from our implementation does not match PyNaCl's result.")

    def test_compare_with_pynacl_shared_secret(self):
        # Generate two random private keys.
        private1 = os.urandom(32)
        private2 = os.urandom(32)
        
        # Create two DiffieHellman objects.
        dh1 = EllipticCurveDiffieHellman(private_key=private1)
        dh2 = EllipticCurveDiffieHellman(private_key=private2)
        
        # Compute shared secret from our implementation.
        shared1 = dh1.generate_shared_secret(dh2.public_key)
        shared2 = dh2.generate_shared_secret(dh1.public_key)
        
        # Compute shared secret using PyNaCl:
        py_shared1 = crypto_scalarmult(private1, dh2.public_key)
        py_shared2 = crypto_scalarmult(private2, dh1.public_key)
        
        # Verify that our implementation and PyNaCl produce the same shared secret.
        self.assertEqual(shared1, py_shared1,
                         msg="Shared secret from our implementation (party 1) does not match PyNaCl's result.")
        self.assertEqual(shared2, py_shared2,
                         msg="Shared secret from our implementation (party 2) does not match PyNaCl's result.")
        self.assertEqual(shared1, shared2,
                         msg="The two computed shared secrets do not match each other.")

    def _test_compare_with_pynacl_box_transformation(self):
        """
        Verify that our shared key (derived via our DH abstraction, which applies
        crypto_core_hsalsa20 to the raw shared secret) matches the key computed by
        PyNaCl's crypto_box_beforenm.
        """
        from nacl.bindings import crypto_box_beforenm
        
        # Generate key pairs using PyNaCl.
        pynacl_a = PrivateKey.generate()
        pynacl_b = PrivateKey.generate()
        
        # Create our Diffie–Hellman instances using the same keys.
        dh_a = EllipticCurveDiffieHellman(private_key=bytes(pynacl_a))
        dh_b = EllipticCurveDiffieHellman(private_key=bytes(pynacl_b))
        
        # Compute the shared key using our implementation.
        # Our generate_shared_secret method clamps, masks, computes the raw shared secret, 
        # then applies crypto_core_hsalsa20 (with sigma) to produce the final key.
        shared_key = dh_a.generate_shared_secret(dh_b.public_key)
        hashed_key = Salsa20.new(shared_key, shared_key[:8])
        # Compute the precomputed key using PyNaCl's crypto_box_beforenm.
        expected_key = crypto_box_beforenm(bytes(pynacl_a), bytes(pynacl_b.public_key))

        self.assertEqual(shared_key, expected_key,
                        msg="Derived shared key does not match PyNaCl's crypto_box_beforenm result.")



if __name__ == '__main__':
    unittest.main()