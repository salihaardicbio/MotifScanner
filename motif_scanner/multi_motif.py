"""
multi_motif.py

Scans a protein database against several motifs at once and
produces one combined table of results.

Motifs may have different lengths. For each database sequence,
every motif is scanned independently (using its own window length)
and its single best-scoring window is kept. Those per-motif best
hits are then combined into one row per protein, so scores stay
comparable across motifs of any length.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

import logging
from pathlib import Path

from motif_scanner.motif_builder import MotifBuilder
from motif_scanner.scanner import Scanner

logger = logging.getLogger(__name__)


def load_motifs_from_directory(motifs_dir, pseudocount=None):
    """
    Builds one MotifProfile per aligned FASTA file found in a
    directory.

    Parameters
    ----------
    motifs_dir
        Directory containing one aligned motif FASTA file per motif.
        The motif's name is taken from the filename (without
        extension). Motifs do not need to be the same length.

    pseudocount
        Optional pseudocount override applied to every motif.

    Returns
    -------
    dict[str, MotifProfile]
        Motif name -> MotifProfile, in filename-sorted order.
    """

    motifs_dir = Path(motifs_dir)

    fasta_files = sorted(
        list(motifs_dir.glob("*.fasta")) + list(motifs_dir.glob("*.fa"))
    )

    if not fasta_files:
        raise ValueError(
            f"No .fasta/.fa files found in '{motifs_dir}'."
        )

    builder = (
        MotifBuilder(pseudocount)
        if pseudocount is not None
        else MotifBuilder()
    )

    profiles = {}

    for fasta_file in fasta_files:

        name = fasta_file.stem

        profiles[name] = builder.build(fasta_file)

    return profiles


class MultiMotifScanner:
    """
    Scans protein sequences against several motif profiles at the
    same time. Motifs may have different lengths - each one is
    scanned independently, and only its single best-scoring window
    per protein is kept.
    """

    def __init__(
        self,
        profiles,
        property_file,
        matrix="BLOSUM80",
        blosum_weight=0.60,
        property_weight=0.20,
        exact_weight=0.20,
    ):
        """
        Parameters
        ----------
        profiles
            dict[str, MotifProfile] - motif name -> profile. Motifs
            may have different lengths.

        property_file, matrix, blosum_weight, property_weight,
        exact_weight
            Passed through to a Scanner built for each motif.
        """

        if not profiles:
            raise ValueError("At least one motif profile is required.")

        self.motif_names = list(profiles.keys())
        self.profiles = profiles

        self.scanners = {
            name: Scanner(
                profile=profile,
                property_file=property_file,
                matrix=matrix,
                blosum_weight=blosum_weight,
                property_weight=property_weight,
                exact_weight=exact_weight,
            )
            for name, profile in profiles.items()
        }

    ##################################################################

    def scan_record(self, record):
        """
        Finds each motif's single best-scoring window in one
        Bio.SeqRecord and combines them into one row.

        Returns
        -------
        dict | None
            None if any motif has no valid window in this sequence
            (e.g. the sequence is shorter than that motif, or
            contains non-standard amino acids) - the whole record is
            skipped in that case, since a fair combined comparison
            needs a result from every motif.
        """

        best_hits = {}

        for name in self.motif_names:

            scanner = self.scanners[name]

            hits = scanner.scan_record(record, minimum_average_score=float("-inf"))

            if not hits:
                logger.warning(
                    "Skipping '%s': no valid window for motif '%s' "
                    "(sequence may be shorter than the motif, or "
                    "contain non-standard amino acids).",
                    record.id,
                    name,
                )
                return None

            best_hits[name] = max(hits, key=lambda hit: hit.average_score)

        row = {
            "sequence_id": record.id,
            "description": record.description,
        }

        multiplied_average = 1.0
        multiplied_total = 1.0

        for name in self.motif_names:

            hit = best_hits[name]

            row[f"{name}_start"] = hit.start
            row[f"{name}_end"] = hit.end
            row[f"{name}_peptide"] = hit.peptide
            row[f"{name}_total"] = hit.total_score
            row[f"{name}_average"] = hit.average_score

            multiplied_average *= hit.average_score
            multiplied_total *= hit.total_score

        row["multiplied_average"] = multiplied_average
        row["multiplied_total"] = multiplied_total

        return row

    ##################################################################

    def scan_records(self, records):
        """
        Scans an iterable of Bio.SeqRecord objects.

        Returns
        -------
        list[dict]
            One row per protein that had a valid window for every
            motif.
        """

        rows = []

        for record in records:

            row = self.scan_record(record)

            if row is not None:
                rows.append(row)

        return rows
