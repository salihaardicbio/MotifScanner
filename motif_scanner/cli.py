"""
cli.py

Command-line interface for MotifScanner.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line argument parser.

    Returns
    -------
    argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(

        prog="MotifScanner",

        description=(
            "Scan protein databases using a motif profile "
            "built from aligned protein sequences."
        )

    )

    ##############################################################
    # Input
    ##############################################################

    input_group = parser.add_argument_group(
        "Input"
    )

    input_group.add_argument(

        "--database",

        required=True,

        type=Path,

        help="Protein FASTA database."

    )

    ##############################################################
    # Motif input
    ##############################################################

    motif_group = parser.add_argument_group(
        "Motif"
    )

    motif = motif_group.add_mutually_exclusive_group(
        required=True
    )

    motif.add_argument(

        "--alignment",

        type=Path,

        help="Aligned motif FASTA."

    )

    motif.add_argument(

        "--profile",

        type=Path,

        help="Saved motif profile (.json)."

    )

    motif.add_argument(

        "--motifs-dir",

        type=Path,

        help=(
            "Directory containing multiple aligned motif FASTA files "
            "(one motif per file, filename used as motif name). Scans "
            "the database against every motif at once (motifs may be "
            "different lengths) and writes a single Excel report "
            "comparing each protein's best match to every motif."
        )

    )

    ##############################################################
    # Output
    ##############################################################

    output_group = parser.add_argument_group(
        "Output"
    )

    output_group.add_argument(

        "--output",

        required=True,

        type=Path,

        help="Output directory."

    )

    output_group.add_argument(

        "--prefix",

        default="motif_scan",

        help="Prefix for output files."

    )

    ##############################################################
    # Scanning options
    ##############################################################

    scan_group = parser.add_argument_group(
        "Scanning"
    )

    scan_group.add_argument(

        "--minimum-score",

        type=float,

        default=0.0,

        help="Minimum average score."

    )

    scan_group.add_argument(

        "--top-hits",

        type=int,

        default=None,

        help="Keep only the top N hits."

    )

    ##############################################################
    # Miscellaneous
    ##############################################################

    parser.add_argument(

        "--verbose",

        action="store_true",

        help="Enable verbose logging."

    )

    parser.add_argument(

        "--version",

        action="version",

        version="MotifScanner 1.0.0"

    )

    return parser


######################################################################


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = build_parser()

    return parser.parse_args()