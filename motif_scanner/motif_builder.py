"""
motif_builder.py

Creates a MotifProfile from an aligned FASTA file.

Pipeline
--------
FASTA
    ↓
Validation
    ↓
Frequency calculation
    ↓
Pseudocount correction
    ↓
Entropy calculation
    ↓
Conservation calculation
    ↓
MotifProfile

Author
------
Saliha Ardıç
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from Bio import SeqIO

from motif_scanner.models import (
    MotifProfile,
    MotifPosition
)

from motif_scanner.constants import DEFAULT_PSEUDOCOUNT

from motif_scanner.pseudocount import (
    PseudocountCalculator
)

from motif_scanner.conservation import (
    calculate_entropy,
    calculate_conservation
)


class MotifBuilder:
    """
    Builds a MotifProfile from an aligned FASTA file.
    """

    def __init__(self, pseudocount: float = DEFAULT_PSEUDOCOUNT):

        self.pseudocount = PseudocountCalculator(pseudocount)

    ##################################################################

    def build(
        self,
        fasta_file
    ):

        fasta_file = Path(fasta_file)

        sequences = [

            str(record.seq).upper()

            for record in SeqIO.parse(
                fasta_file,
                "fasta"
            )

        ]

        ##############################################################

        if len(sequences) == 0:

            raise ValueError(
                "No sequences found."
            )

        ##############################################################

        length = len(sequences[0])

        for sequence in sequences:

            if len(sequence) != length:

                raise ValueError(
                    "Sequences must already be aligned."
                )

        ##############################################################

        profile = MotifProfile(

            length=length,

            sequence_count=len(sequences),

            positions=[]

        )

        ##############################################################

        for position in range(length):

            residues = [

                sequence[position]

                for sequence in sequences

            ]

            ##########################################################

            counts = Counter(residues)

            # Raw observed frequencies (used only for entropy/conservation)
            raw_frequencies = {
                aa: count / len(residues)
                for aa, count in counts.items()
            }

            # Pseudocount frequencies (used for scoring)
            frequencies = self.pseudocount.calculate(
                residues
            )

            ##########################################################

            consensus = max(

                counts,

                key=counts.get

            )

            ##########################################################

            entropy = calculate_entropy(
                raw_frequencies
            )
            
            conservation = calculate_conservation(
                entropy
            )

            ##########################################################

            motif_position = MotifPosition(

                index=position + 1,

                consensus=consensus,

                observed_residues=set(
                    residues
                ),

                counts=dict(
                    counts
                ),

                frequencies=frequencies,

                entropy=entropy,

                conservation_weight=conservation

            )

            ##########################################################

            profile.positions.append(
                motif_position
            )

        ##############################################################

        return profile