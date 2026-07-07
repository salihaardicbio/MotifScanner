import unittest

from motif_scanner.scanner import Scanner


class TestScanner(unittest.TestCase):

    def test_clean_sequence(self):

        sequence = "acdefg"

        cleaned = Scanner.clean_sequence(
            sequence
        )

        self.assertEqual(

            cleaned,

            "ACDEFG"

        )

    def test_validate_sequence(self):

        self.assertTrue(

            Scanner.validate_sequence(

                "ACDEFGHIKLMNPQRSTVWY"

            )

        )

        self.assertFalse(

            Scanner.validate_sequence(

                "ACDZ"

            )

        )


if __name__ == "__main__":

    unittest.main()