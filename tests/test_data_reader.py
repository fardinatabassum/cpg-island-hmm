import unittest
from src.data_reader import read_fasta, get_header
import os

class TestDataReader(unittest.TestCase):
    def setUp(self):
        # Create a small, messy temporary FASTA file
        self.test_file = "test_temp.fa"
        with open(self.test_file, "w") as f:
            f.write(">chr22 [source=GRCh38] Homo sapiens\n")
            f.write("ATGC\n")
            f.write("atgc\n") # Testing case sensitivity
            f.write("NNNN\n") # Testing gaps

    def tearDown(self):
        # Clean up the file after testing
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_header_cleaning(self):
        """Verifies that only the ID (chr22) is extracted, not the metadata."""
        header = get_header(self.test_file)
        self.assertEqual(header, "chr22")

    def test_sequence_continuity(self):
        """Verifies that newlines are stripped and case is unified."""
        sequence = read_fasta(self.test_file)
        # Expected: ATGC + atgc + NNNN = ATGCATGCNNNN (12bp)
        self.assertEqual(len(sequence), 12)
        self.assertEqual(sequence[:8].upper(), "ATGCATGC")

if __name__ == "__main__":
    unittest.main()