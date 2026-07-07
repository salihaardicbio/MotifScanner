"""
scanner.py

Slides a motif profile across protein sequences and scores
every window.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

import logging

from motif_scanner.constants import AMINO_ACIDS, WINDOW_STEP
from motif_scanner.models import WindowResult
from motif_scanner.scorer import PositionScorer

logger = logging.getLogger(__name__)

STANDARD_AMINO_ACIDS = set(AMINO_ACIDS)


class Scanner:
    """
    Scans protein sequences against a MotifProfile using a
    PositionScorer.
    """

    def __init__(
        self,
        profile,
        property_file,
        matrix="BLOSUM80",
        blosum_weight=0.60,
        property_weight=0.20,
        exact_weight=0.20,
    ):
        """
        Parameters
        ----------
        profile
            MotifProfile object built by MotifBuilder or loaded
            via ProfileIO.

        property_file
            Path to the amino acid physicochemical properties CSV.

        matrix
            Substitution matrix name (e.g. "BLOSUM80").

        blosum_weight, property_weight, exact_weight
            Component weights used by the PositionScorer.
        """

        self.profile = profile
        self.motif_length = profile.length

        self.scorer = PositionScorer(
            property_file,
            matrix=matrix,
            blosum_weight=blosum_weight,
            property_weight=property_weight,
            exact_weight=exact_weight,
        )

    ##################################################################

    @staticmethod
    def clean_sequence(sequence: str) -> str:
        """
        Normalizes a raw sequence string.

        Strips whitespace and converts to uppercase.
        """

        return sequence.strip().upper()

    ##################################################################

    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """
        Returns True only if every residue is one of the 20
        standard amino acids.
        """

        if not sequence:
            return False

        return all(
            residue in STANDARD_AMINO_ACIDS
            for residue in sequence
        )

    ##################################################################

    def _scan_sequence(
        self,
        sequence_id,
        description,
        sequence,
        minimum_average_score=0.0,
    ):
        """
        Slides the motif window across a single cleaned, validated
        sequence and scores every window.

        Returns
        -------
        list[WindowResult]
        """

        results = []

        sequence_length = len(sequence)

        if sequence_length < self.motif_length:
            return results

        last_start = sequence_length - self.motif_length

        for start in range(0, last_start + 1, WINDOW_STEP):

            end = start + self.motif_length

            peptide = sequence[start:end]

            position_scores = []
            total_score = 0.0

            for offset, motif_position in enumerate(self.profile.positions):

                query_residue = peptide[offset]

                position_score = self.scorer.score(
                    query_residue,
                    motif_position,
                )

                position_scores.append(position_score)
                total_score += position_score.final

            average_score = total_score / self.motif_length

            if average_score >= minimum_average_score:

                results.append(
                    WindowResult(
                        sequence_id=sequence_id,
                        description=description,
                        start=start + 1,
                        end=end,
                        peptide=peptide,
                        total_score=total_score,
                        average_score=average_score,
                        position_scores=position_scores,
                    )
                )

        return results

    ##################################################################

    def scan_record(
        self,
        record,
        minimum_average_score=0.0
    ):
        """
        Scan a single Bio.SeqRecord object.

        Parameters
        ----------
        record
            Bio.SeqRecord object.

        minimum_average_score : float

        Returns
        -------
        list[WindowResult]
        """

        sequence = self.clean_sequence(
            str(record.seq)
        )

        if not self.validate_sequence(
            sequence
        ):
            logger.warning(
                "Skipping '%s' because it contains "
                "non-standard amino acids.",
                record.id
            )
            return []

        return self._scan_sequence(
            sequence_id=record.id,
            description=record.description,
            sequence=sequence,
            minimum_average_score=minimum_average_score
        )

    ##################################################################

    def scan_records(
        self,
        records,
        minimum_average_score=0.0
    ):
        """
        Scan an iterable of SeqRecord objects.

        Parameters
        ----------
        records
            Iterable of Bio.SeqRecord objects.

        minimum_average_score : float

        Returns
        -------
        list[WindowResult]
        """

        results = []

        for record in records:

            results.extend(

                self.scan_record(

                    record,

                    minimum_average_score

                )

            )

        return results

    ##################################################################

    def summary(self, results, proteins_scanned, windows_scored, elapsed_time):
        """
        Builds a summary dictionary describing a completed scan.

        Parameters
        ----------
        results
            Final (filtered) list of WindowResult objects (hits).

        proteins_scanned
            Number of database sequences that were scanned.

        windows_scored
            Total number of sliding windows that were scored,
            before filtering by minimum_average_score.

        elapsed_time
            Wall-clock time in seconds.

        Returns
        -------
        dict
        """

        best_score = max(
            (result.average_score for result in results),
            default=0.0,
        )

        return {
            "proteins": proteins_scanned,
            "windows": windows_scored,
            "hits": len(results),
            "best_score": best_score,
            "elapsed_time": elapsed_time,
        }

    ##################################################################

    def __repr__(self):
        """
        String representation of the Scanner.
        """

        return (

            f"{self.__class__.__name__}("

            f"motif_length={self.motif_length}, "

            f"sequence_count={self.profile.sequence_count})"

        )


##########################################################################
# End of file
##########################################################################
