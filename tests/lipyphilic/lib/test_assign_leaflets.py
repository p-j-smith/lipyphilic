
import pytest
import numpy as np
import MDAnalysis

from numpy.testing import assert_array_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT
from lipyphilic.lib.assign_leaflets import AssignLeaflets
 
 
class TestAssignLeaflets:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C'
    }
    
    @pytest.fixture(scope='class')
    def leaflets(self, universe):
        leaflets = AssignLeaflets(universe, **self.kwargs)
        leaflets.run()
        return leaflets

    def test_assign_leaflets(self, leaflets):
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'leaflets_present': [-1, 1]  # all lipids should be assigned to the lower (-1) or upper (1) leaflet
        }
    
        assert leaflets.leaflets.shape == (reference['n_residues'], reference['n_frames'])
        assert (np.unique(leaflets.leaflets) == reference['leaflets_present']).all()
    
        reference = {
            'assigned': np.array([[1]] * 50 + [[-1]] * 50)  # first 50 residues are in the upper leaflet
        }

        assert_array_equal(leaflets.leaflets, reference['assigned'])


class TestAssignLeafletsExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    def test_Exceptions(self, universe):
        
        match = "midplane_sel contains atoms that are not present in molecules selected "
        with pytest.raises(ValueError, match=match):
            AssignLeaflets(
                universe=universe,
                lipid_sel="name L",
                midplane_sel="name C",
                midplane_cutoff=10
            )
