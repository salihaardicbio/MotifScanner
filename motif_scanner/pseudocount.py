"""
pseudocount.py

Applies Bayesian pseudocount correction
to amino acid frequencies.
"""

from collections import Counter

from motif_scanner.constants import (
    AMINO_ACIDS,
    BACKGROUND_FREQUENCIES,
    DEFAULT_PSEUDOCOUNT,
)


class PseudocountCalculator:

    def __init__(self,
                 pseudocount=DEFAULT_PSEUDOCOUNT):

        self.lambda_ = pseudocount

    #######################################################################

    def calculate(self, residues):

        """
        residues

        Example

        ["I","I","I","G"]
        """

        counts = Counter(residues)

        N = len(residues)

        frequencies = {}

        denominator = N + self.lambda_

        for aa in AMINO_ACIDS:

            observed = counts.get(aa, 0)

            background = BACKGROUND_FREQUENCIES[aa]

            frequencies[aa] = (
                observed +
                self.lambda_ * background
            ) / denominator

        return frequencies