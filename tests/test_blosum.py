import unittest

from motif_scanner.blosum import BlosumMatrix


class TestBlosum(unittest.TestCase):

    def setUp(self):

        self.matrix = BlosumMatrix()

    def test_identity(self):

        score = self.matrix.score(
            "A",
            "A"
        )

        self.assertGreaterEqual(
            score,
            0.0
        )

        self.assertLessEqual(
            score,
            1.0
        )

    def test_substitution(self):

        score = self.matrix.score(
            "A",
            "V"
        )

        self.assertGreaterEqual(
            score,
            0.0
        )

        self.assertLessEqual(
            score,
            1.0
        )


if __name__ == "__main__":

    unittest.main()