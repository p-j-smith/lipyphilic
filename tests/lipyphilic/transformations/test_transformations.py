import MDAnalysis
import numpy as np
import pytest

from lipyphilic._simple_systems.simple_systems import (
    HEX_LAT_SPLIT_Z,
    TRICLINIC,
)
from lipyphilic.transformations import (
    CentreMembrane,
)


class TestCenterMembrane:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_SPLIT_Z)

    reference = {
        "upper_leaflet": np.arange(100),
        "lower_leaflet": np.arange(100, 200),
        "broken_bilayer_height": 92,
        "broken_bilayer_midpoint": 0,
        "correct_bilayer_height": 8,
        "bilayer_midpoint": 50,
    }

    def test_no_center(self, universe):
        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])

        assert bilayer_height == self.reference["broken_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]

    def test_center_frame(self, universe):
        membrane = universe.select_atoms("name L C")
        ts = universe.trajectory[0]
        CentreMembrane(ag=membrane, shift=5)(ts)

        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])

        assert bilayer_height == self.reference["correct_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]

    def test_center_trajectory(self):
        universe = MDAnalysis.Universe(HEX_LAT_SPLIT_Z)

        membrane = universe.select_atoms("name L C")
        universe.trajectory.add_transformations(CentreMembrane(ag=membrane))

        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])

        assert bilayer_height == self.reference["correct_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]

    def test_exceptions(self):
        universe_triclinic = MDAnalysis.Universe(TRICLINIC)

        match = "CentreMembrane requires an orthorhombic box"
        with pytest.raises(ValueError, match=match):
            universe_triclinic.trajectory.add_transformations(
                CentreMembrane(ag=universe_triclinic.atoms),
            )
