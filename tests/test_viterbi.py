import unittest
import math
import numpy as np
from src.viterbi_hmm import viterbi_decode

class TestViterbi(unittest.TestCase):
    def setUp(self):
        self.states = [0, 1]
        self.start_p = [math.log(0.5), math.log(0.5)]
        self.trans_p = {0: {0: math.log(0.95), 1: math.log(0.05)}, 1: {0: math.log(0.05), 1: math.log(0.95)}}
        self.emiss_p = {0: {'A': math.log(0.25), 'C': math.log(0.25), 'G': math.log(0.25), 'T': math.log(0.25)},
                        1: {'A': math.log(0.01), 'C': math.log(0.49), 'G': math.log(0.49), 'T': math.log(0.01)}}

    def test_n_gap_alignment(self):
        """Verifies coordinates are preserved across N-gaps."""
        seq = "G" * 20 + "N" * 10 + "G" * 20
        path = viterbi_decode(seq, self.states, self.start_p, self.trans_p, self.emiss_p)
        self.assertEqual(len(path), 50) # Should be 50 base pairs

    def test_background_bias(self):
        """Verifies AT-rich sequences stay in background state."""
        seq = "ATATATATAT" * 5 # 50bp of AT
        path = viterbi_decode(seq, self.states, self.start_p, self.trans_p, self.emiss_p)
        self.assertTrue(all(s == 0 for s in path)) # Should all be background (0)

if __name__ == "__main__":
    unittest.main()
