"""
MotifScanner

A motif scanning package using BLOSUM80, physicochemical
properties, and conservation-weighted scoring.
"""

__version__ = "1.0.0"

from .scanner import Scanner
from .motif_builder import MotifBuilder
from .multi_motif import MultiMotifScanner, load_motifs_from_directory
from .profile_io import ProfileIO
from .report import ReportWriter

__all__ = [
    "Scanner",
    "MotifBuilder",
    "MultiMotifScanner",
    "load_motifs_from_directory",
    "ProfileIO",
    "ReportWriter",
]