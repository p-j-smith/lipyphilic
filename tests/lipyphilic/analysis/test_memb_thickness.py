import MDAnalysis
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from lipyphilic._simple_systems.simple_systems import HEX_LAT, HEX_LAT_BUMP, TRICLINIC
from lipyphilic.analysis.memb_thickness import MembThickness


class TestMembThickness:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    leaflets = np.ones((100, 1))
    leaflets[50:] = -1

    kwargs = {
        "lipid_sel": "name L C",
        "leaflets": leaflets,
    }

    def test_memb_thickness(self, universe):
        memb_thickness = MembThickness(universe, **self.kwargs)
        memb_thickness.run()

        reference = {
            "n_frames": 1,
            "thickness": [20],
        }

        assert memb_thickness.results.memb_thickness.shape == (reference["n_frames"],)
        assert_array_equal(memb_thickness.results.memb_thickness, reference["thickness"])


class TestMembThicknessUndulating:
    @staticmethod
    @pytest.fixture
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP)

    leaflets = np.ones((100, 1))
    leaflets[50:] = -1

    kwargs = {
        "lipid_sel": "name L C",
        "leaflets": leaflets,
    }

    reference = {
        "thickness": [20],
        "thickness_grid": np.array(
            [
                [
                    [30.0, 30.0, 30.0, 30.0],
                    [20.0, 20.0, 20.0, 20.0],
                    [20.0, 20.0, 20.0, 20.0],
                    [20.0, 20.0, 20.0, 20.0],
                ],
            ],
        ),  # 1 frame with 4 by 4 membrane grid
    }

    def test_nbins4(self, universe):
        memb_thickness = MembThickness(universe, n_bins=4, **self.kwargs)
        memb_thickness.run()

        assert_array_equal(memb_thickness.results.memb_thickness, self.reference["thickness"])

    def test_nbins200_no_interpolation(self, universe):
        memb_thickness = MembThickness(universe, n_bins=200, **self.kwargs)
        memb_thickness.run()

        reference = {
            "thickness": [np.nan],
        }

        assert_array_equal(memb_thickness.results.memb_thickness, reference["thickness"])

    def test_nbins200_with_interpolation(self, universe):
        memb_thickness = MembThickness(universe, n_bins=200, interpolate=True, **self.kwargs)
        memb_thickness.run()

        reference = {
            "thickness": [20],
        }

        assert_array_equal(memb_thickness.results.memb_thickness, reference["thickness"])

    def test_return_surface(self, universe):
        # Shift the first 20 atoms by 10 Å in z so that the thickness is non-uniform
        # We need to move the universe into memory otherwise the translation will be lost
        # when iterating over frames
        universe.transfer_to_memory()
        universe.atoms[:20].translate((0, 0, 10))
        memb_thickness = MembThickness(universe, n_bins=4, return_surface=True, **self.kwargs)
        memb_thickness.run()

        assert_array_equal(memb_thickness.memb_thickness_grid, self.reference["thickness_grid"])


class TestMembThicknessExceptions:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    def test_Exceptions(self, universe):
        match = "'leaflets' must either be a 1D array containing non-changing "
        with pytest.raises(ValueError, match=match):
            MembThickness(
                universe=universe,
                lipid_sel="name L C",
                leaflets=None,
            )

        with pytest.raises(ValueError, match=match):
            MembThickness(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[[], []], [[], []]]),  # cannot pass a 3D array
            )

        match = (
            "The shape of 'leaflets' must be \\(n_residues,\\), but 'lipid_sel' "
            "generates an AtomGroup with 100 residues"
            " and 'leaflets' has shape \\(99, 1\\)."
        )
        with pytest.raises(ValueError, match=match):
            MembThickness(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[1]] * 50 + [[-1]] * 49),  # one residue too few
            )

        match = "The frames to analyse must be identical to those used in assigning lipids to leaflets."
        with pytest.raises(ValueError, match=match):
            areas = MembThickness(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[1, 1]] * 50 + [[-1, -1]] * 50),  # leaflets has two frames, apl one
            )
            areas.run()

        universe_triclinic = MDAnalysis.Universe(TRICLINIC)
        match = "MembThickness requires an orthorhombic box"
        with pytest.raises(ValueError, match=match):
            MembThickness(
                universe=universe_triclinic,
                lipid_sel="name C",
                leaflets=[0, 0],
            )
