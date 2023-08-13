import MDAnalysis
import numpy as np
from numpy.testing import assert_array_equal
import pytest

import lipyphilic as lpp
from lipyphilic._simple_systems.simple_systems import (
    HEX_LAT,
    HEX_LAT_BUMP,
    HEX_LAT_BUMP_MID_ATOM,
    HEX_LAT_BUMP_MID_MOL,
    TRICLINIC,
)


class TestAssignLeaflets:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        "lipid_sel": "name L C",
    }

    @pytest.fixture(scope="class")
    def leaflets(self, universe):
        leaflets = lpp.AssignLeaflets(universe, **self.kwargs)
        leaflets.run()
        return leaflets

    def test_assign_leaflets(self, leaflets):
        reference = {
            "n_residues": 100,
            "n_frames": 1,
            "leaflets_present": [
                -1,
                1,
            ],  # all lipids should be assigned to the lower (-1) or upper (1) leaflet
        }

        assert leaflets.results.leaflets.shape == (reference["n_residues"], reference["n_frames"])
        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])

        # first 50 residues are in the upper leaflet (1)
        # final 50 residues are in the lower leaflet (-1)
        reference = {
            "assigned": np.array([[1]] * 50 + [[-1]] * 50),
        }

        assert_array_equal(leaflets.results.leaflets, reference["assigned"])

    def test_filter_leaflets(self, leaflets):
        reference = {
            "n_residues": 50,
            "n_frames": 1,
            "leaflets_present": [-1, 1],
        }

        filtered_leaflets = leaflets.filter_leaflets(lipid_sel="name C", start=None, stop=None, step=None)

        assert filtered_leaflets.shape == (reference["n_residues"], reference["n_frames"])
        assert_array_equal(np.unique(filtered_leaflets), reference["leaflets_present"])

        # first 25 residues are in the upper leaflet (1)
        # final 25 residues are in the lower leaflet (-1)
        reference = {
            "assigned": np.array([[1]] * 25 + [[-1]] * 25),
        }

        assert_array_equal(filtered_leaflets, reference["assigned"])


class TestAssignLeafletsExceptions:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    def test_Exceptions(self, universe):
        match = "midplane_sel is 'None' and midplane_cutoff "
        with pytest.raises(ValueError, match=match):
            lpp.AssignLeaflets(
                universe=universe,
                lipid_sel="name L",
                midplane_sel=None,
                midplane_cutoff=5,
            )

        match = "midplane_sel is 'name C' and midplane_cutoff "
        with pytest.raises(ValueError, match=match):
            lpp.AssignLeaflets(
                universe=universe,
                lipid_sel="name L C",
                midplane_sel="name C",
                midplane_cutoff=None,
            )

        match = "To assign molecules to the midplane, midplane_cutoff must"
        with pytest.raises(ValueError, match=match):
            lpp.AssignLeaflets(
                universe=universe,
                lipid_sel="name L C",
                midplane_sel="name C",
                midplane_cutoff=-10,
            )

        match = "midplane_sel contains atoms that are not present in molecules selected "
        with pytest.raises(ValueError, match=match):
            lpp.AssignLeaflets(
                universe=universe,
                lipid_sel="name L",
                midplane_sel="name C",
                midplane_cutoff=10,
            )

        universe_triclinic = MDAnalysis.Universe(TRICLINIC)
        match = "AssignLeaflets requires an orthorhombic box. Please use the on-the-fly"
        with pytest.raises(ValueError, match=match):
            lpp.AssignLeaflets(
                universe=universe_triclinic,
                lipid_sel="name C",
            )


class TestAssignLeafletsUndulating:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP)

    kwargs = {
        "lipid_sel": "name L C",
        "midplane_sel": "name C",
        "midplane_cutoff": 6.5,
    }

    def test_nbins1(self, universe):
        leaflets = lpp.AssignLeaflets(universe, n_bins=1, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 0, 1],
            "midplane_resnames": ["CHOL"] * 6,  # list or residues incorrectly identified as midplane
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )
        assert "LIPID" not in universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames

    def test_nbins4(self, universe):
        leaflets = lpp.AssignLeaflets(universe, n_bins=4, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 1],  # now (correctly) no midplane residues should be found
            "midplane_resnames": [],  # list or residues incorrectly identified as midplane
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )


class TestAssignLeafletsUndulatingMidplaneMol:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP_MID_MOL)

    kwargs = {
        "lipid_sel": "name L C",
        "midplane_sel": "name C",
        "midplane_cutoff": 6.5,
        "n_bins": 4,
    }

    @pytest.fixture(scope="class")
    def leaflets(self, universe):
        leaflets = lpp.AssignLeaflets(universe, **self.kwargs)
        leaflets.run()
        return leaflets

    def test_nbins4_midplane(self, universe):
        leaflets = lpp.AssignLeaflets(universe, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 0, 1],
            "midplane_resnames": ["CHOL"],  # list of residues resnames (correctly) identified as midplane
            "midplane_resids": [78],  # resid of midplane molecules
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resids,
            reference["midplane_resids"],
        )


class TestAssignLeafletsUndulatingMidplaneAtom:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP_MID_ATOM)

    kwargs = {
        "lipid_sel": "name L C",
        "midplane_sel": "name C",
        "midplane_cutoff": 6.5,
        "n_bins": 4,
    }

    @pytest.fixture(scope="class")
    def leaflets(self, universe):
        leaflets = lpp.AssignLeaflets(universe, **self.kwargs)
        leaflets.run()
        return leaflets

    def test_nbins4_midplane(self, universe):
        leaflets = lpp.AssignLeaflets(universe, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 1],
            "midplane_resnames": [],  # only one of two atoms in CHOL78 is in the midplane, so the molecule is not assigned to the midplane
            "midplane_resids": [],
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resids,
            reference["midplane_resids"],
        )


class TestAssignCurvedLeafletsUndulating:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP)

    kwargs = {
        "lipid_sel": "name L C",
        "lf_cutoff": 12,
    }

    def test_no_midplane(self, universe):
        leaflets = lpp.AssignCurvedLeaflets(universe, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 1],
            "midplane_resnames": [],
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )


class TestAssignCurvedLeafletsUndulatingMidplaneMol:
    @staticmethod
    @pytest.fixture(scope="class")
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP_MID_MOL)

    kwargs = {
        "lipid_sel": "name L C",
        "lf_cutoff": 12,
        "midplane_sel": "name C",
        "midplane_cutoff": 12,
    }

    def test_midplane(self, universe):
        leaflets = lpp.AssignCurvedLeaflets(universe, **self.kwargs)
        leaflets.run()

        reference = {
            "leaflets_present": [-1, 0, 1],
            "midplane_resnames": ["CHOL"],
            "midplane_resids": [78],  # resid of midplane molecules
        }

        assert_array_equal(np.unique(leaflets.results.leaflets), reference["leaflets_present"])
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resnames,
            reference["midplane_resnames"],
        )
        assert_array_equal(
            universe.residues[leaflets.results.leaflets[:, 0] == 0].resids,
            reference["midplane_resids"],
        )
