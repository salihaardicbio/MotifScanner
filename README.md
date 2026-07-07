# MotifScanner

A protein motif scanning tool that combines BLOSUM substitution scores,
physicochemical property similarity, and conservation-weighted scoring
to find matches to a motif profile in a protein database.

Author: Saliha Ardıç

## How it works

1. **Build a motif profile** from an aligned set of example sequences
   (a FASTA file where every sequence has the same length). For each
   column of the alignment, MotifScanner computes amino acid
   frequencies (with pseudocount correction), Shannon entropy, and a
   conservation weight.
2. **Scan a protein database** by sliding a window the length of the
   motif across every sequence. Each window is scored against the
   motif profile using a weighted combination of:
   - **BLOSUM substitution score** (default: BLOSUM80)
   - **Physicochemical property similarity** (hydrophobicity, volume,
     polarity, charge, flexibility, aromaticity, helix/sheet
     propensity)
   - **Exact match bonus**
   - A **conservation weight**, which down-weights positions that are
     highly variable in the original alignment.
3. **Reports** are written as CSV, JSON, and plain-text summaries.

## Installation

```bash
pip install -r requirements.txt
pip install -e .
```

This installs the `motif-scanner` command-line tool.

## Project Structure

```text
MotifScanner/
├── motif_scanner/
│   ├── models.py                  # Data classes (MotifProfile, WindowResult, ...)
│   ├── motif_builder.py           # Builds a MotifProfile from an aligned FASTA
│   ├── scorer.py                  # Scores one residue against one motif position
│   ├── scanner.py                 # Slides the motif across a database and scores windows
│   ├── profile_io.py              # Save/load motif profiles as JSON
│   ├── report.py                  # Writes CSV/JSON/text reports
│   ├── cli.py                     # Command-line argument parsing
│   ├── main.py                    # Entry point (ties everything together)
│   ├── blosum.py                  # Normalized substitution matrix scoring
│   ├── conservation.py            # Shannon entropy and conservation weights
│   ├── pseudocount.py             # Pseudocount-corrected amino acid frequencies
│   ├── property_similarity.py     # Physicochemical similarity between amino acids
│   ├── amino_acid_properties.py   # Loads and normalizes aa_properties.csv
│   ├── multi_motif.py             # Scans a database against several motifs at once (any lengths)
│   └── aa_properties.csv          # Amino acid physicochemical descriptor table
├── config.yaml                    # Default scoring weights, matrix, logging, etc.
├── requirements.txt
├── pyproject.toml
└── README.md
```

## Usage

### Build a profile from an alignment and scan a database

```bash
motif-scanner \
  --database examples/database.fasta \
  --alignment examples/motif.fasta \
  --output results/
```

### Scan using a previously saved profile

```bash
motif-scanner \
  --database examples/database.fasta \
  --profile results/motif_profile.json \
  --output results/
```

### Scan a database against multiple motifs at once (Excel report)

If you have several motifs and want to compare how a database scores
against all of them side by side, put one aligned FASTA file per
motif into a folder (the filename becomes the motif's name in the
report) and use `--motifs-dir` instead of `--alignment`/`--profile`:

```bash
motif-scanner \
  --database examples/database.fasta \
  --motifs-dir motifs/ \
  --output results/
```

**Motifs can be different lengths.** Each motif is scanned
independently against every protein, and only that motif's single
best-scoring window is kept. Those best hits are then combined into
one row per protein. If a protein is shorter than one of the motifs
(or otherwise has no valid window for it), that protein is skipped
for the whole comparison, with a warning logged, since a fair
combined score needs a result from every motif.

This produces a single Excel file,
`<prefix>_multi_motif_results.xlsx`, with one row per protein:
`Sequence_ID`, `Description`, then for every motif a
`<name> - Start`, `<name> - End`, `<name> - Peptide`,
`<name> - Total_Score`, and `<name> - Average_Score` column group
(that motif's own best window in this protein), and finally two
combined columns — `Multiplied_Average_Score` and
`Multiplied_Total_Score` — which are the product of every motif's
best average (respectively total) score for that protein. A protein
that scores well against every motif will have a much higher
multiplied score than one that only matches a single motif, since a
single low score pulls the whole product down. Rows are sorted by
`Multiplied_Average_Score`, highest first.

### Options

| Flag | Description | Default |
|---|---|---|
| `--database` | Protein FASTA database to scan (required) | — |
| `--alignment` | Aligned motif FASTA (mutually exclusive with `--profile`/`--motifs-dir`) | — |
| `--profile` | Saved motif profile JSON (mutually exclusive with `--alignment`/`--motifs-dir`) | — |
| `--motifs-dir` | Folder of aligned motif FASTA files, one per motif — motifs can be different lengths (mutually exclusive with `--alignment`/`--profile`) — produces an Excel comparison report instead of the usual CSV/JSON/text reports | — |
| `--output` | Output directory for reports (required) | — |
| `--prefix` | Prefix for output filenames | `motif_scan` |
| `--minimum-score` | Minimum average score for a hit to be kept (in `--motifs-dir` mode, a protein's row is kept if **at least one** motif's best average score meets this threshold) | `0.0` |
| `--top-hits` | Keep only the top N hits | all hits |
| `--verbose` | Enable verbose logging and a progress bar | off |
| `--version` | Print the version and exit | — |

`--alignment`, `--profile`, and `--motifs-dir` are mutually exclusive:
use `--alignment` to build a fresh profile, `--profile` to reuse one
saved earlier, or `--motifs-dir` to compare several motifs at once.

## Configuration

`config.yaml` sets defaults that aren't exposed as CLI flags:

```yaml
scoring:
  matrix: BLOSUM80
  weights:
    blosum: 0.60
    property: 0.20
    exact: 0.20

motif:
  pseudocount: 1.0

logging:
  level: INFO
```

## Output files

Each scan produces, in the output directory:

- `<prefix>_results.csv` — one row per hit
- `<prefix>_results.json` — full results, including per-position scores
- `<prefix>_summary.txt` — scan summary (proteins scanned, windows scored, hits, best score, elapsed time)
- `<prefix>_position_details.csv` — per-position score breakdown for every hit
- `<prefix>_top_hits.txt` — human-readable top hits report

When using `--motifs-dir`, instead of the above you get a single
`<prefix>_multi_motif_results.xlsx` comparing every motif (see
"Scan a database against multiple motifs at once" above).

## Running tests

```bash
pytest tests/
```

## License

MIT
