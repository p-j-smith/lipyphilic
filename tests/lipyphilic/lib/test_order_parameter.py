
import pytest
import numpy as np
import MDAnalysis

from numpy.testing._private.utils import assert_array_almost_equal, assert_array_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT
from lipyphilic.lib.order_parameter import SCC
 
 
class TestSCC:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'tail_sel': 'name L C',
    }
    
    def test_SCC(self, universe):
        
        scc = SCC(universe, **self.kwargs)
        scc.run()
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }
        
        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])
        
    def test_SCC_normals(self, universe):
        
        normals = np.ones((100, 1))
        
        scc = SCC(universe, **self.kwargs, normals=normals)
        scc.run()
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }
        
        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])
        
    def test_SCC_normals_3D(self, universe):
            
        normals = np.zeros((100, 1, 3))
        normals[:, :, 2] = 1.0
        
        scc = SCC(universe, **self.kwargs, normals=normals)
        scc.run()
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }
        
        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])


class TestSCCWeightedAverage:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)
    
    @pytest.fixture(scope='class')
    def sn1_scc(self, universe):
        sn1_scc = SCC(universe, "name L")
        sn1_scc.run()
        return sn1_scc
    
    @pytest.fixture(scope='class')
    def sn2_scc(self, universe):
        sn2_scc = SCC(universe, "name C")
        sn2_scc.run()
        return sn2_scc
    
    def test_SCC_weighted_average(self, sn1_scc):
        
        scc = SCC.weighted_average(sn1_scc, sn1_scc)
        
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'scc': np.full((50, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }

        assert scc.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc, reference['scc'])
        
    def test_SCC_weighted_average_different_tails(self, sn1_scc, sn2_scc):
        
        scc = SCC.weighted_average(sn1_scc, sn2_scc)
        
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }

        assert scc.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc, reference['scc'])
        
    def test_SCC_weighted_average_indices(self, sn1_scc, sn2_scc):
        
        scc, indices = SCC.weighted_average(sn1_scc, sn2_scc, return_indices=True)
        
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5),  # all bonds are perpendicular to the z-axis
            'indices': np.arange(100)
        }

        assert scc.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc, reference['scc'])
        assert_array_equal(indices, reference['indices'])


class TestSCCExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'tail_sel': 'name L C',
    }

    def test_Exceptions(self, universe):
            
        match = "'normals' must either be a 2D array containing leaflet ids "
        with pytest.raises(ValueError, match=match):
            SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones(100)
            )
        
        match = "The shape of 'normals' must be \\(n_residues,\\)"
        with pytest.raises(ValueError, match=match):
            SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones((99, 1))
            )
            
        match = "The frames to analyse must be identical to those used "
        with pytest.raises(ValueError, match=match):
            scc = SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones((100, 2))
            )
            scc.run()
            

class TestSCCWeightedAverageExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'tail_sel': 'name L C',
    }
    
    @pytest.fixture(scope='class')
    def sn1_scc(self, universe):
        sn1_scc = SCC(universe, **self.kwargs)
        sn1_scc.run()
        return sn1_scc
    
    def test_Exceptions(self, universe, sn1_scc):
            
        match = "sn1_scc and sn2_scc must have been run with the same frames"
        with pytest.raises(ValueError, match=match):
            sn2_scc = SCC(
                universe=universe,
                **self.kwargs,
            )
            sn2_scc.run(stop=0)
            SCC.weighted_average(sn1_scc, sn2_scc)
            
        with pytest.raises(ValueError, match=match):
            sn2_scc = SCC(
                universe=universe,
                **self.kwargs,
            )
            sn2_scc.run(stop=1)
            sn2_scc.frames = np.array([10])
            SCC.weighted_average(sn1_scc, sn2_scc)
