import tempfile
import unittest

from motif_scanner.profile_io import ProfileIO


class TestProfileIO(unittest.TestCase):

    def test_validate(self):

        with tempfile.NamedTemporaryFile(
            suffix=".json"
        ) as f:

            self.assertFalse(

                ProfileIO.validate(
                    f.name
                )

            )


if __name__ == "__main__":

    unittest.main()