"""
constants.py

Global constants used throughout the MotifScanner project.

Author : Saliha Ardıç
Project : Property-Augmented Position Specific Motif Scanner (PA-PSSM)

This module stores:

- Standard amino acids
- Background amino acid frequencies
- Default scoring weights
- Default pseudocount value
- Supported substitution matrices
- AAindex properties used in the project
"""

from typing import List, Dict

###############################################################################
# Amino Acids
###############################################################################

AMINO_ACIDS: List[str] = [
    "A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V"
]

###############################################################################
# Background frequencies
#
# Robinson & Robinson (1991)
# These are commonly used as amino acid background frequencies.
###############################################################################

BACKGROUND_FREQUENCIES: Dict[str, float] = {

    "A":0.07805,
    "R":0.05129,
    "N":0.04487,
    "D":0.05364,
    "C":0.01925,
    "Q":0.04264,
    "E":0.06295,
    "G":0.07377,
    "H":0.02199,
    "I":0.05142,
    "L":0.09019,
    "K":0.05744,
    "M":0.02243,
    "F":0.03856,
    "P":0.05203,
    "S":0.07120,
    "T":0.05841,
    "W":0.01330,
    "Y":0.03216,
    "V":0.06441
}

###############################################################################
# Default pseudocount
###############################################################################

DEFAULT_PSEUDOCOUNT: float = 1.0

###############################################################################
# Default scoring weights
###############################################################################

DEFAULT_SCORING_WEIGHTS = {

    "blosum":0.70,

    "property":0.20,

    "exact_match":0.10
}

###############################################################################
# Supported substitution matrices
###############################################################################

SUPPORTED_MATRICES = [

    "BLOSUM80",

    "BLOSUM62",

    "PAM250"
]

###############################################################################
# AAindex descriptors
#
# These descriptors were selected because they describe different
# physicochemical aspects of amino acids while minimizing redundancy.
###############################################################################

AA_PROPERTIES = [

    "hydrophobicity",

    "side_chain_volume",

    "polarity",

    "net_charge",

    "flexibility",

    "aromaticity",

    "helix_propensity",

    "sheet_propensity"
]

###############################################################################
# Default configuration
###############################################################################

DEFAULT_MATRIX = "BLOSUM80"

DEFAULT_OUTPUT_DIRECTORY = "results"

DEFAULT_PROFILE_NAME = "motif_profile.json"

DEFAULT_REPORT_DIRECTORY = "reports"

###############################################################################
# Window scanning
###############################################################################

WINDOW_STEP = 1

###############################################################################
# Exact match score
###############################################################################

EXACT_MATCH_SCORE = 1.0

NO_MATCH_SCORE = 0.0