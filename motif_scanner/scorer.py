"""
scorer.py

Calculates motif position scores.

Author
------
Saliha Ardıç
"""

from __future__ import annotations

from motif_scanner.property_similarity import PropertySimilarity
from motif_scanner.blosum import BlosumMatrix
from motif_scanner.models import PositionScore


class PositionScorer:
    """
    Scores a single residue against one motif position.
    """

    def __init__(
        self,
        property_file,
        matrix="BLOSUM80",
        blosum_weight=0.60,
        property_weight=0.20,
        exact_weight=0.20
    ):

        self.property_similarity = PropertySimilarity(
            property_file
        )

        self.blosum = BlosumMatrix(matrix)

        self.w_blosum = blosum_weight
        self.w_property = property_weight
        self.w_exact = exact_weight

    ##############################################################

    def score(
        self,
        query_residue,
        motif_position
    ):
        """
        Score one query residue against one motif position.

        Parameters
        ----------
        query_residue : str
            Amino acid from the database sequence.

        motif_position : MotifPosition
            One position from the motif profile.

        Returns
        -------
        PositionScore
        """

        ##########################################################
        # BLOSUM component
        ##########################################################

        blosum_score = 0.0

        for residue, frequency in motif_position.frequencies.items():

            blosum_score += (

                frequency *

                self.blosum.score(
                    query_residue,
                    residue
                )

            )


        ##########################################################
        # Physicochemical similarity
        ##########################################################

        property_score = 0.0

        for residue, frequency in motif_position.frequencies.items():

            property_score += (

                frequency *

                self.property_similarity.similarity(
                    query_residue,
                    residue
                )

            )

        ##########################################################
        # Exact match
        ##########################################################

        if query_residue in motif_position.observed_residues:

            exact_score = 1.0

        else:

            exact_score = 0.0

        ##########################################################
        # Conservation
        ##########################################################

        conservation_score = motif_position.conservation_weight

        ##########################################################
        # Final score
        ##########################################################

        base_score = (

            self.w_blosum * blosum_score +

            self.w_property * property_score +

            self.w_exact * exact_score

        )

        final_score = (

            base_score *

            conservation_score

        )

        ##########################################################


        return PositionScore(

            position=motif_position.index,

            query_residue=query_residue,

            blosum=blosum_score,

            property=property_score,

            exact=exact_score,

            conservation=conservation_score,

            final=final_score

        )