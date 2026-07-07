"""
models.py

Core data models used throughout MotifScanner.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set


# ============================================================
# Position within a motif profile
# ============================================================

@dataclass
class MotifPosition:
    """
    Represents one column of the motif alignment.
    """

    index: int

    consensus: str

    observed_residues: Set[str]

    counts: Dict[str, int]

    frequencies: Dict[str, float]

    entropy: float

    conservation_weight: float


# ============================================================
# Complete motif profile
# ============================================================

@dataclass
class MotifProfile:
    """
    Complete motif profile.
    """

    length: int

    sequence_count: int

    positions: List[MotifPosition]


# ============================================================
# Score of a single motif position
# ============================================================

@dataclass
class PositionScore:
    """
    Score of one motif position.
    """

    position: int

    query_residue: str

    blosum: float

    property: float

    exact: float

    conservation: float

    final: float


# ============================================================
# One motif hit
# ============================================================

@dataclass
class WindowResult:
    """
    One sliding-window hit.
    """

    sequence_id: str

    description: str

    start: int

    end: int

    peptide: str

    total_score: float

    average_score: float

    position_scores: List[PositionScore] = field(
        default_factory=list
    )