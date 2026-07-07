"""
main.py

Entry point for MotifScanner.

Ties together the CLI, motif building/loading, scanning,
and report generation.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import yaml
from Bio import SeqIO
from tqdm import tqdm

from motif_scanner.cli import parse_arguments
from motif_scanner.constants import DEFAULT_MATRIX
from motif_scanner.motif_builder import MotifBuilder
from motif_scanner.multi_motif import MultiMotifScanner, load_motifs_from_directory
from motif_scanner.profile_io import ProfileIO
from motif_scanner.report import ReportWriter
from motif_scanner.scanner import Scanner

PACKAGE_DIRECTORY = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_DIRECTORY.parent

DEFAULT_CONFIG_FILE = PROJECT_ROOT / "config.yaml"
DEFAULT_PROPERTY_FILE = PACKAGE_DIRECTORY / "aa_properties.csv"


def load_config(config_file: Path = DEFAULT_CONFIG_FILE) -> dict:
    """
    Loads config.yaml if it exists.

    Returns an empty dict (falling back to built-in defaults)
    if the file is missing.
    """

    if not config_file.exists():
        return {}

    with open(config_file, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def configure_logging(verbose: bool, level_name: str) -> logging.Logger:
    """
    Configures root logging for the CLI run.
    """

    level = logging.DEBUG if verbose else getattr(
        logging, str(level_name).upper(), logging.INFO
    )

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    return logging.getLogger("motif_scanner")


def run_multi_motif_scan(args, config, logger) -> int:
    """
    Scans the database against every motif found in --motifs-dir and
    writes a single Excel comparison report.

    Motifs may have different lengths: each one is scanned
    independently and only its best-scoring window per protein is
    kept, then combined into one row per protein.
    """

    scoring_config = config.get("scoring", {})
    weights = scoring_config.get("weights", {})
    matrix_name = scoring_config.get("matrix", DEFAULT_MATRIX)

    motif_config = config.get("motif", {})
    pseudocount = motif_config.get("pseudocount")

    logger.info("Loading motifs from directory: %s", args.motifs_dir)

    profiles = load_motifs_from_directory(args.motifs_dir, pseudocount=pseudocount)

    logger.info(
        "Loaded %d motif(s): %s",
        len(profiles),
        ", ".join(profiles.keys()),
    )

    scanner = MultiMotifScanner(
        profiles=profiles,
        property_file=str(DEFAULT_PROPERTY_FILE),
        matrix=matrix_name,
        blosum_weight=weights.get("blosum", 0.60),
        property_weight=weights.get("property", 0.20),
        exact_weight=weights.get("exact", 0.20),
    )

    logger.info("Scanning database: %s", args.database)

    records = list(SeqIO.parse(args.database, "fasta"))

    all_rows = []

    total_records = len(records)

    progress = tqdm(records, desc="Scanning database", unit="seq")

    for index, record in enumerate(progress, start=1):

        progress.set_postfix(remaining=total_records - index)

        row = scanner.scan_record(record)

        if row is not None:
            all_rows.append(row)

    all_rows = [
        row for row in all_rows
        if any(
            row[f"{name}_average"] >= args.minimum_score
            for name in scanner.motif_names
        )
    ]

    all_rows.sort(key=lambda row: row["multiplied_average"], reverse=True)

    if args.top_hits is not None:
        all_rows = all_rows[: args.top_hits]

    args.output.mkdir(parents=True, exist_ok=True)

    output_file = args.output / f"{args.prefix}_multi_motif_results.xlsx"

    ReportWriter.write_multi_motif_excel(
        all_rows,
        motif_names=scanner.motif_names,
        output_file=output_file,
    )

    logger.info("Wrote %d row(s) to: %s", len(all_rows), output_file)

    print()
    print("=" * 60)
    print("Generated Report File")
    print("=" * 60)
    print(f"{'excel':20s} : {output_file}")
    print("=" * 60)

    return 0


def main() -> int:
    """
    MotifScanner command-line entry point.
    """

    args = parse_arguments()
    config = load_config()

    logger = configure_logging(
        args.verbose,
        config.get("logging", {}).get("level", "INFO"),
    )

    if args.motifs_dir is not None:
        return run_multi_motif_scan(args, config, logger)

    scoring_config = config.get("scoring", {})
    weights = scoring_config.get("weights", {})
    matrix_name = scoring_config.get("matrix", DEFAULT_MATRIX)

    motif_config = config.get("motif", {})
    pseudocount = motif_config.get("pseudocount")

    # --------------------------------------------------------------
    # Build or load the motif profile.
    # --------------------------------------------------------------

    if args.alignment is not None:

        logger.info("Building motif profile from alignment: %s", args.alignment)

        builder = (
            MotifBuilder(pseudocount=pseudocount)
            if pseudocount is not None
            else MotifBuilder()
        )

        profile = builder.build(args.alignment)

    else:

        logger.info("Loading motif profile from: %s", args.profile)

        profile = ProfileIO.load(args.profile)

    logger.info(
        "Motif length: %d, sequences used: %d",
        profile.length,
        profile.sequence_count,
    )

    # --------------------------------------------------------------
    # Scan the database.
    # --------------------------------------------------------------

    scanner = Scanner(
        profile=profile,
        property_file=str(DEFAULT_PROPERTY_FILE),
        matrix=matrix_name,
        blosum_weight=weights.get("blosum", 0.60),
        property_weight=weights.get("property", 0.20),
        exact_weight=weights.get("exact", 0.20),
    )

    logger.info("Scanning database: %s", args.database)

    records = list(SeqIO.parse(args.database, "fasta"))

    start_time = time.time()

    all_results = []
    windows_scored = 0

    total_records = len(records)

    progress = tqdm(records, desc="Scanning database", unit="seq")

    for index, record in enumerate(progress, start=1):

        progress.set_postfix(remaining=total_records - index)

        sequence = Scanner.clean_sequence(str(record.seq))

        if Scanner.validate_sequence(sequence) and len(sequence) >= profile.length:
            windows_scored += len(sequence) - profile.length + 1

        all_results.extend(
            scanner.scan_record(
                record,
                minimum_average_score=args.minimum_score,
            )
        )

    elapsed_time = time.time() - start_time

    all_results.sort(key=lambda result: result.average_score, reverse=True)

    if args.top_hits is not None:
        all_results = all_results[: args.top_hits]

    summary = scanner.summary(
        all_results,
        proteins_scanned=len(records),
        windows_scored=windows_scored,
        elapsed_time=elapsed_time,
    )

    # --------------------------------------------------------------
    # Write reports.
    # --------------------------------------------------------------

    files = ReportWriter.write_all(
        all_results,
        summary,
        args.output,
        prefix=args.prefix,
    )

    ReportWriter.print_report_locations(files)

    return 0


if __name__ == "__main__":
    sys.exit(main())
