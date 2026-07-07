"""
profile_io.py

Save and load MotifProfile objects as JSON.

This module is responsible only for serialization and
deserialization. It does NOT build motif profiles or perform
scoring.

Author:
    Saliha Ardıç
"""

from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime

from .models import MotifProfile, MotifPosition


SOFTWARE_NAME = "MotifScanner"
SOFTWARE_VERSION = "1.0.0"


class ProfileIO:
    """
    Read and write MotifProfile objects.

    Notes
    -----
    The JSON file stores only biological information and metadata.
    Computational lookup tables are intentionally NOT stored so they
    can be regenerated in future versions.
    """

    ##################################################################
    @staticmethod
    def save(profile: MotifProfile, output_file: str | Path) -> None:
        """
        Save a MotifProfile to a JSON file.

        Parameters
        ----------
        profile
            MotifProfile object.
        output_file
            Destination JSON filename.
        """

        output_file = Path(output_file)

        data = {
            "software": {
                "name": SOFTWARE_NAME,
                "version": SOFTWARE_VERSION
            },
            "created": datetime.now().isoformat(timespec="seconds"),
            "motif_length": profile.length,
            "sequence_count": profile.sequence_count,
            "positions": []
        }

        for position in profile.positions:

            data["positions"].append({

                "index": position.index,

                "consensus": position.consensus,

                "observed_residues":
                    sorted(list(position.observed_residues)),

                "counts":
                    position.counts,

                "frequencies":
                    position.frequencies,

                "entropy":
                    position.entropy,

                "conservation_weight":
                    position.conservation_weight

            })

        with open(output_file, "w", encoding="utf-8") as handle:

            json.dump(
                data,
                handle,
                indent=4
            )

    ##################################################################
    @staticmethod
    def load(profile_file: str | Path) -> MotifProfile:
        """
        Load a MotifProfile from JSON.

        Parameters
        ----------
        profile_file
            JSON profile.

        Returns
        -------
        MotifProfile
        """

        profile_file = Path(profile_file)

        with open(profile_file, "r", encoding="utf-8") as handle:

            data = json.load(handle)

        positions = []

        for p in data["positions"]:

            position = MotifPosition(

                index=p["index"],

                consensus=p["consensus"],

                observed_residues=set(
                    p["observed_residues"]
                ),

                counts=p["counts"],

                frequencies=p["frequencies"],

                entropy=p["entropy"],

                conservation_weight=p["conservation_weight"]

            )

            positions.append(position)

        profile = MotifProfile(

            length=data["motif_length"],

            sequence_count=data["sequence_count"],

            positions=positions

        )

        return profile

    ##################################################################
    @staticmethod
    def info(profile_file: str | Path) -> dict:
        """
        Read only profile metadata without loading
        the complete profile.

        Useful for quickly checking a profile.

        Returns
        -------
        dict
        """

        profile_file = Path(profile_file)

        with open(profile_file, "r", encoding="utf-8") as handle:

            data = json.load(handle)

        return {

            "software":
                data.get("software", {}),

            "created":
                data.get("created"),

            "motif_length":
                data.get("motif_length"),

            "sequence_count":
                data.get("sequence_count")

        }

    ##################################################################
    @staticmethod
    def validate(profile_file: str | Path) -> bool:
        """
        Basic validation of a profile JSON.

        Returns
        -------
        bool
        """

        try:

            info = ProfileIO.info(profile_file)

            required = [

                "software",

                "motif_length",

                "sequence_count"

            ]

            for key in required:

                if key not in info:

                    return False

            return True

        except Exception:

            return False