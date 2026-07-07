import unittest

from motif_scanner.conservation import Conservation


class TestConservation(unittest.TestCase):

    def test_full_conservation(self):

        frequencies = {

            "A": 1.0

        }

        score = Conservation.weight(
            frequencies
        )

        self.assertAlmostEqual(
            score,
            1.0,
            places=6
        )

    def test_partial_conservation(self):

        frequencies = {

            "A": 0.5,

            "G": 0.5

        }

        score = Conservation.weight(
            frequencies
        )

        self.assertGreater(
            score,
            0.0
        )

        self.assertLess(
            score,
            1.0
        )


if __name__ == "__main__":

    unittest.main()