import random
import unittest

from tonellishanks import tonellishanks

from tonelli import tonelli


class TestTonelliShanks(unittest.TestCase):
    def _test_tonelli_shanks_against_naive(self) -> None:
        p = 2**255 - 19
        # We only need this to work for 2^255 - 19, lol
        # But should work WLOG (trust me, bro)

        for _ in range(2, 1_000):
            n = random.randint(1, 10_000_000)
            ts = tonellishanks(n, p)
            t = ts**2 % p if ts is not None else None
            res = tonelli(n, p)
            t2 = res**2 % p if res is not None else None
            self.assertEqual(t, t2, f"n: {n}, p: {p}")
            self.assertEqual(ts, res, f"n: {n}, p: {p}")
