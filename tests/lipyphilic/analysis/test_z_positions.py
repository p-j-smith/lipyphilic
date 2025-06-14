import MDAnalysis
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from lipyphilic._simple_systems.simple_systems import HEX_LAT, HEX_LAT_BUMP, TRICLINIC
from lipyphilic.analysis.z_positions import ZPositions


class TestZPositions:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        "lipid_sel": "name L C",
        "height_sel": "name C",
        "n_bins": 1,
    }

    @pytest.fixture(scope="class")
    def z_positions(self, universe):
        z_positions = ZPositions(universe, **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
        reference = {
            "n_residues": 50,
            "n_frames": 1,
            "z_positions": np.full((50, 1), fill_value=10),
        }

        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference["z_positions"][25:] = -10

        assert z_positions.results.z_positions.shape == (reference["n_residues"], reference["n_frames"])
        assert_array_equal(z_positions.results.z_positions, reference["z_positions"])

    def test_exceptions(self):
        universe_triclinic = MDAnalysis.Universe(TRICLINIC)
        match = "ZPositions requires an orthorhombic box"
        with pytest.raises(ValueError, match=match):
            ZPositions(
                universe=universe_triclinic,
                lipid_sel="name C",
                height_sel="name C C",
            )


class TestZPositionsOneAtom:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        "lipid_sel": "name L C",
        "n_bins": 1,
    }

    @pytest.fixture(scope="class")
    def z_positions(self, universe):
        height_indices = " ".join(universe.atoms[2::4].indices.astype(str))  # one atom per cholesterol
        z_positions = ZPositions(universe, height_sel=f"index {height_indices}", **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
        reference = {
            "n_residues": 50,
            "n_frames": 1,
            "z_positions": np.full((50, 1), fill_value=10),
        }

        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference["z_positions"][25:] = -10

        assert z_positions.results.z_positions.shape == (reference["n_residues"], reference["n_frames"])
        assert_array_equal(z_positions.results.z_positions, reference["z_positions"])


class TestZPositionsUndulating:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP)

    kwargs = {
        "lipid_sel": "name L C",
        "height_sel": "name C",
        "n_bins": 10,  # with fewer bins the z positions will no longer be +/- 10 Angstrom
    }

    @pytest.fixture(scope="class")
    def z_positions(self, universe):
        z_positions = ZPositions(universe, **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
        reference = {
            "n_residues": 50,
            "n_frames": 1,
            "z_positions": np.full((50, 1), fill_value=10),
        }

        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference["z_positions"][25:] = -10

        assert z_positions.results.z_positions.shape == (reference["n_residues"], reference["n_frames"])
        assert_array_equal(z_positions.results.z_positions, reference["z_positions"])

    def test_return_memb_midpoint(self, universe):
        reference = {
            "memb_midpoint": np.full((1, 10, 10), fill_value=50),
        }
        reference["memb_midpoint"][0, 3:7, 2:8] = 55
        z_positions = ZPositions(universe, return_midpoint=True, **self.kwargs)
        z_positions.run()
        assert z_positions.memb_midpoint.shape == reference["memb_midpoint"].shape
        assert_array_equal(z_positions.memb_midpoint, reference["memb_midpoint"])
