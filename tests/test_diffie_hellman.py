import os
import binascii
import unittest
from parameterized import parameterized

from nacl.bindings import crypto_scalarmult

from x25519.elliptic_curve_diffie_hellman import EllipticCurveDiffieHellman
from x25519.group_law import X25519CurveGroupLaw
from x25519.montgomery_ladder import MontgomeryLadderMKTutorial, MontgomeryLadderOptimized, MontgomeryLadderRFC7748

class TestDiffieHellmanVectors(unittest.TestCase):
    def setUp(self):
        self.base_point = b'\x09' + (b'\x00' * 31)

    @parameterized.expand([
        ("MontgomeryLadderRFC7748", MontgomeryLadderRFC7748()),
        ("MontgomeryLadderMKTutorial", MontgomeryLadderMKTutorial()),
        ("MontgomeryLadderOptimized", MontgomeryLadderOptimized()),
        ("GroupLaw", X25519CurveGroupLaw()),
    ])
    def test_alice_public_key(self, name, curve):
        # Test vector from RFC 7748
        # Alice's private key, a:
        alice_private_hex = "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a"
        # Expected Alice's public key, X25519(a, 9):
        expected_alice_public_hex = "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a"
        
        alice_private = binascii.unhexlify(alice_private_hex)
        expected_alice_public = binascii.unhexlify(expected_alice_public_hex)
        
        alice = EllipticCurveDiffieHellman(private_key=alice_private, curve=curve)
        self.assertEqual(alice.public_key, expected_alice_public,
                         msg=f"Alice public key does not match the test vector for {name}.")

    @parameterized.expand([
        ("MontgomeryLadderRFC7748", MontgomeryLadderRFC7748()),
        ("MontgomeryLadderMKTutorial", MontgomeryLadderMKTutorial()),
        ("MontgomeryLadderOptimized", MontgomeryLadderOptimized()),
        ("GroupLaw", X25519CurveGroupLaw()),
    ])
    def test_bob_public_key(self, name, curve):
        # Bob's private key, b:
        bob_private_hex = "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb"
        # Expected Bob's public key, X25519(b, 9):
        expected_bob_public_hex = "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f"
        
        bob_private = binascii.unhexlify(bob_private_hex)
        expected_bob_public = binascii.unhexlify(expected_bob_public_hex)
        
        bob = EllipticCurveDiffieHellman(private_key=bob_private, curve=curve)
        self.assertEqual(bob.public_key, expected_bob_public,
                         msg=f"Bob public key does not match the test vector for {name}.")

    @parameterized.expand([
        ("MontgomeryLadderRFC7748", MontgomeryLadderRFC7748()),
        ("MontgomeryLadderMKTutorial", MontgomeryLadderMKTutorial()),
        ("MontgomeryLadderOptimized", MontgomeryLadderOptimized()),
        ("GroupLaw", X25519CurveGroupLaw()),
    ])
    def test_shared_secret(self, name, curve):
        # Using the same test vectors as above:
        alice_private_hex = "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a"
        bob_private_hex = "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb"
        # Expected shared secret, K:
        expected_shared_hex = "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742"
        
        alice_private = binascii.unhexlify(alice_private_hex)
        bob_private = binascii.unhexlify(bob_private_hex)
        expected_shared = binascii.unhexlify(expected_shared_hex)
        
        alice = EllipticCurveDiffieHellman(private_key=alice_private, curve=curve)
        bob = EllipticCurveDiffieHellman(private_key=bob_private, curve=curve)
        
        # Alice computes the shared secret using Bob's public key.
        alice_shared = alice.generate_shared_secret(bob.public_key)
        # Bob computes the shared secret using Alice's public key.
        bob_shared = bob.generate_shared_secret(alice.public_key)
        
        self.assertEqual(alice_shared, expected_shared,
                         msg=f"Alice's computed shared secret does not match the test vector for {name}.")
        self.assertEqual(bob_shared, expected_shared,
                         msg=f"Bob's computed shared secret does not match the test vector for {name}.")
        self.assertEqual(alice_shared, bob_shared,
                         msg=f"Alice and Bob's shared secrets do not match for {name}.")

    @parameterized.expand([
        ("MontgomeryLadderRFC7748", MontgomeryLadderRFC7748()),
        ("MontgomeryLadderMKTutorial", MontgomeryLadderMKTutorial()),
        ("MontgomeryLadderOptimized", MontgomeryLadderOptimized()),
        ("GroupLaw", X25519CurveGroupLaw()),
    ])
    def test_compare_with_pynacl_public_key(self, name, curve):
        # Generate a random private key.
        private = os.urandom(32)
        
        # Compute public key using our DiffieHellman abstraction.
        my_dh = EllipticCurveDiffieHellman(private_key=private, curve=curve)
        my_public = my_dh.public_key
        
        # Compute public key using PyNaCl's crypto_scalarmult.
        py_public = crypto_scalarmult(private, self.base_point)
        
        self.assertEqual(my_public, py_public,
                         msg=f"Public key from our implementation does not match PyNaCl's result for {name}.")

    @parameterized.expand([
        ("MontgomeryLadderRFC7748", MontgomeryLadderRFC7748()),
        ("MontgomeryLadderMKTutorial", MontgomeryLadderMKTutorial()),
        ("MontgomeryLadderOptimized", MontgomeryLadderOptimized()),
        ("GroupLaw", X25519CurveGroupLaw()),
    ])
    def test_compare_with_pynacl_shared_secret(self, name, curve):
        # Generate two random private keys.
        private1 = os.urandom(32)
        private2 = os.urandom(32)
        
        # Create two DiffieHellman objects.
        dh1 = EllipticCurveDiffieHellman(private_key=private1, curve=curve)
        dh2 = EllipticCurveDiffieHellman(private_key=private2, curve=curve)
        
        # Compute shared secret from our implementation.
        shared1 = dh1.generate_shared_secret(dh2.public_key)
        shared2 = dh2.generate_shared_secret(dh1.public_key)
        
        # Compute shared secret using PyNaCl:
        py_shared1 = crypto_scalarmult(private1, dh2.public_key)
        py_shared2 = crypto_scalarmult(private2, dh1.public_key)
        
        self.assertEqual(shared1, py_shared1,
                         msg=f"Shared secret from our implementation (party 1) does not match PyNaCl's result for {name}.")
        self.assertEqual(shared2, py_shared2,
                         msg=f"Shared secret from our implementation (party 2) does not match PyNaCl's result for {name}.")
        self.assertEqual(shared1, shared2,
                         msg=f"The two computed shared secrets do not match each other for {name}.")

    """
    def _test_compare_with_pynacl_box_transformation(self):
     
        #Verify that our shared key (derived via our DH abstraction, which applies
        #crypto_core_hsalsa20 to the raw shared secret) matches the key computed by
        #PyNaCl's crypto_box_beforenm.
        from nacl.bindings import crypto_box_beforenm
        
        # Generate key pairs using PyNaCl.
        pynacl_a = PrivateKey.generate()
        pynacl_b = PrivateKey.generate()
        
        # Create our Diffie-Hellman instances using the same keys.
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
    """


if __name__ == '__main__':
    unittest.main()