"""
blosum.py

Loads substitution matrices and returns normalized scores.

Normalization
-------------
Scores are normalized per query residue (row-wise), not against the
matrix's global minimum and maximum. For a given query residue
`aa1`, the score is scaled relative to the best and worst raw
substitution score `aa1` can achieve against any standard amino
acid:

score_norm(aa1, aa2) =
    (raw_score(aa1, aa2) - row_min[aa1]) / (row_max[aa1] - row_min[aa1])

This keeps every score in [0, 1] while using the full range for each
residue's own realistic best/worst match, instead of compressing
most substitutions into whatever narrow band is left over by the
matrix's single most extreme pair.

Author
------
Saliha Ardıç
"""

from Bio.Align import substitution_matrices

from motif_scanner.constants import AMINO_ACIDS


class BlosumMatrix:

    def __init__(self, matrix_name="BLOSUM80"):

        self.matrix = substitution_matrices.load(matrix_name)

        self.row_min, self.row_max = self._calculate_row_ranges()

    ##################################################################

    def _calculate_row_ranges(self):
        """
        Precomputes, for every standard amino acid, the minimum and
        maximum raw substitution score it can achieve against any
        other standard amino acid.
        """

        row_min = {}
        row_max = {}

        for aa in AMINO_ACIDS:

            scores = [
                self.raw_score(aa, other)
                for other in AMINO_ACIDS
            ]

            row_min[aa] = min(scores)
            row_max[aa] = max(scores)

        return row_min, row_max

    ##################################################################

    def raw_score(
        self,
        aa1,
        aa2
    ):

        try:

            return self.matrix[(aa1, aa2)]

        except KeyError:

            return self.matrix[(aa2, aa1)]

    ##################################################################

    def score(
        self,
        aa1,
        aa2
    ):
        """
        Returns a per-residue (row-wise) normalized score.

        aa1 is treated as the query residue: the score is scaled
        relative to aa1's own best and worst possible substitution
        score against any standard amino acid.

        Output

        0 ≤ score ≤ 1
        """

        raw = self.raw_score(
            aa1,
            aa2
        )

        minimum = self.row_min[aa1]
        maximum = self.row_max[aa1]

        if maximum == minimum:
            return 1.0

        normalized = (raw - minimum) / (maximum - minimum)

        return normalized
