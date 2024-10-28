import unittest
from src import needleman_wunsch, smith_waterman

class TestAlignmentAlgorithms(unittest.TestCase):
    def setUp(self):
        # Ortak parametreler
        self.seq1 = "ATAT"
        self.seq2 = "TATA"
        self.substitution_matrix_path = "substitution_matrix.csv"  # CSV dosyasının yolu
        self.gap_penalty = -2
        self.max_alignments = 2
        self.output_file = "output.txt"

    def test_needleman_wunsch(self):
        # Global alignment testi
        score, alignments = needleman_wunsch(self.seq1, self.seq2, self.substitution_matrix_path, self.gap_penalty, self.max_alignments, self.output_file)
        self.assertIsNotNone(score, "Score should not be None")
        self.assertGreaterEqual(len(alignments), 1, "At least one alignment should be returned")
        self.assertTrue(all("-" in alignment or alignment.replace("-", "") in self.seq1 + self.seq2 for alignment in alignments), "Alignments should contain valid characters")

    def test_smith_waterman(self):
        # Local alignment testi
        score, alignments = smith_waterman(self.seq1, self.seq2, self.substitution_matrix_path, self.gap_penalty, self.max_alignments, self.output_file)
        self.assertIsNotNone(score, "Score should not be None")
        self.assertGreaterEqual(len(alignments), 1, "At least one alignment should be returned")
        self.assertTrue(all("-" in alignment or alignment.replace("-", "") in self.seq1 + self.seq2 for alignment in alignments), "Alignments should contain valid characters")

    def test_output_file(self):
        # Dosya çıktısını kontrol etme
        needleman_wunsch(self.seq1, self.seq2, self.substitution_matrix_path, self.gap_penalty, self.max_alignments, self.output_file)
        with open(self.output_file, "r") as file:
            output_content = file.read()
        self.assertIn("Global alignment no.", output_content, "Output file should contain alignment headers")
        self.assertIn("Score:", output_content, "Output file should contain scores for each alignment")

if __name__ == "__main__":
    unittest.main()
