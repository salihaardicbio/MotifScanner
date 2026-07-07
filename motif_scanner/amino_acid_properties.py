"""
amino_acid_properties.py

Loads amino acid physicochemical descriptors and
creates normalized descriptors used during scoring.
"""

from pathlib import Path
from typing import List

import pandas as pd


class AminoAcidProperties:

    def __init__(self, csv_file: str):

        self.csv_file = Path(csv_file)

        self.raw_table = self._load_table()

        self.normalized_table = self._normalize()

    ###########################################################################

    def _load_table(self):

        if not self.csv_file.exists():

            raise FileNotFoundError(
                f"{self.csv_file} does not exist."
            )

        df = pd.read_csv(self.csv_file)

        df = df.set_index("AA")

        return df

    ###########################################################################

    def _normalize(self):

        """
        Z-score normalization.

        Every property becomes

        mean = 0

        std = 1
        """

        df = self.raw_table.copy()

        for column in df.columns:

            mean = df[column].mean()

            std = df[column].std()

            df[column] = (df[column] - mean) / std

        return df

    ###########################################################################

    def get_vector(self, amino_acid: str):

        """
        Normalized descriptor vector.
        """

        return self.normalized_table.loc[amino_acid]

    ###########################################################################

    def get_raw_vector(self, amino_acid: str):

        """
        Original descriptor vector.

        Used only for reporting.
        """

        return self.raw_table.loc[amino_acid]

    ###########################################################################

    def property_names(self) -> List[str]:

        return list(self.raw_table.columns)

    ###########################################################################

    def raw(self):

        return self.raw_table.copy()

    ###########################################################################

    def normalized(self):

        return self.normalized_table.copy()