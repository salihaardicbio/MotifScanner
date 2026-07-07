"""
property_similarity.py

Calculates physicochemical similarity between amino acids.

The comparison is performed independently for each
physicochemical descriptor.

The final similarity is the average of all descriptor
similarities.

Author:
    Saliha Ardıç
"""

from typing import Dict

from motif_scanner.amino_acid_properties import (
    AminoAcidProperties
)


class PropertySimilarity:

    def __init__(self,
                 property_file):

        self.properties = AminoAcidProperties(
            property_file
        )

        self.max_ranges = self._calculate_ranges()

    ##################################################################

    def _calculate_ranges(self):

        """
        Calculates the maximum range of every descriptor.

        Used for normalization.
        """

        table = self.properties.normalized()

        ranges = {}

        for column in table.columns:

            column_range = (

                table[column].max()

                -

                table[column].min()

            )

            # Guard against division by zero if a descriptor is
            # constant across all amino acids.
            ranges[column] = column_range if column_range > 0 else 1.0

        return ranges

    ##################################################################

    def similarity(
        self,
        aa1,
        aa2
    ):

        """
        Calculates similarity between two amino acids.

        Returns

        value between

        0

        and

        1
        """

        vector1 = self.properties.get_vector(aa1)

        vector2 = self.properties.get_vector(aa2)

        similarities = []

        ##############################################################

        for property_name in self.properties.property_names():

            difference = abs(

                vector1[property_name]

                -

                vector2[property_name]

            )

            maximum = self.max_ranges[property_name]

            similarity = 1 - (

                difference

                /

                maximum

            )

            similarities.append(similarity)

        ##############################################################

        return sum(similarities) / len(similarities)