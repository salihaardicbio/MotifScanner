"""
conservation.py

Calculates Shannon entropy and conservation weights
for motif positions.

Author:
    Saliha Ardıç

Project:
    Property-Augmented Position Specific Motif Scanner
"""

import math
from typing import Dict


###############################################################################

def calculate_entropy(
    frequencies: Dict[str, float]
) -> float:
    """
    Calculates Shannon entropy from amino acid frequencies.

    Parameters
    ----------
    frequencies : dict

        Amino acid frequencies.

    Returns
    -------
    float

        Shannon entropy.
    """

    entropy = 0.0

    for probability in frequencies.values():

        if probability > 0:

            entropy -= probability * math.log2(probability)

    return entropy


###############################################################################

def calculate_conservation(
    entropy: float
) -> float:
    """
    Converts entropy into a conservation weight.

    Weight ranges from

    0 = completely variable

    1 = completely conserved
    """

    max_entropy = math.log2(20)

    conservation = 1 - (entropy / max_entropy)

    return conservation


###############################################################################

def calculate_entropy_and_conservation(
    frequencies: Dict[str, float]
):
    """
    Convenience function.

    Returns

    entropy,
    conservation
    """

    entropy = calculate_entropy(frequencies)

    conservation = calculate_conservation(entropy)

    return entropy, conservation


###############################################################################

class Conservation:
    """
    Class-based convenience wrapper around the module-level
    entropy/conservation functions.
    """

    @staticmethod
    def entropy(frequencies: Dict[str, float]) -> float:
        return calculate_entropy(frequencies)

    @staticmethod
    def weight(frequencies: Dict[str, float]) -> float:
        """
        Computes the conservation weight directly from a set of
        amino acid frequencies.
        """

        entropy = calculate_entropy(frequencies)

        return calculate_conservation(entropy)