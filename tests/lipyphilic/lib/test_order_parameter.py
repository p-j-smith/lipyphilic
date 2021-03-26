
import pytest
import numpy as np
import matplotlib
import MDAnalysis

from numpy.testing._private.utils import assert_array_almost_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT, ONE_CHOL, ONE_CHOL_TRAJ
from lipyphilic.lib.order_parameter import SCC
from lipyphilic.lib.plotting import ProjectionPlot
 
matplotlib.use("Agg")
 
 
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
        
    def test_SCC_normals_3D(self, universe):
        
        # The Scc should be the same whether the normal is the positive or negative z-axis
        normals = np.zeros((100, 1, 3))
        normals[:, :50, 2] = 1.0
        normals[:, 50:, 2] = -1.0
        
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

        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])
        
    def test_SCC_weighted_average_different_tails(self, sn1_scc, sn2_scc):
        
        scc = SCC.weighted_average(sn1_scc, sn2_scc)
        
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5)  # all bonds are perpendicular to the z-axis
        }

        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])
        
    def test_SCC_weighted_average_indices(self, sn1_scc, sn2_scc):
        
        scc = SCC.weighted_average(sn1_scc, sn2_scc)
        
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'scc': np.full((100, 1), fill_value=-0.5),  # all bonds are perpendicular to the z-axis
        }

        assert scc.SCC.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(scc.SCC, reference['scc'])


class TestSCCExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'tail_sel': 'name L C',
    }

    def test_Exceptions(self, universe):
            
        match = "'normals' must be a 3D array containing local membrane normals of each lipi at each frame."
        with pytest.raises(ValueError, match=match):
            SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones(100)  # Wrong dimension
            )
        
        match = "The shape of 'normals' must be \\(n_residues, n_frames, 3\\)"
        with pytest.raises(ValueError, match=match):
            SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones((99, 1, 3))  # Too few molecules
            )
            
        match = "The frames to analyse must be identical to those used "
        with pytest.raises(ValueError, match=match):
            scc = SCC(
                universe=universe,
                **self.kwargs,
                normals=np.ones((100, 2, 3))  # Too many frames
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


class TestSCCProjectSCC:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    @pytest.fixture(scope="class")
    def scc(self, universe):
        scc = SCC(
            universe=universe,
            tail_sel="not name ROH"
        )
        scc.run()
        
        return scc
    
    def test_project_SCC(self, scc):
        
        scc_projection = scc.project_SCC()
        
        assert isinstance(scc_projection, ProjectionPlot)
        assert_array_almost_equal(scc_projection.values, scc.SCC.mean())

    def test_filter_by(self, scc):
        
        scc_projection = scc.project_SCC(filter_by=[True])
        
        assert_array_almost_equal(scc_projection.values, scc.SCC.mean())
        
    def test_filter_by_2D(self, scc):
        
        scc_projection = scc.project_SCC(filter_by=np.full((1, 25), fill_value=True))
        
        assert_array_almost_equal(scc_projection.values, scc.SCC.mean())
        
    def test_filter_by_exception(self, scc):
        
        match = "The shape of `filter_by` must either be \(n_lipids, n_frames\) or \(n_lipids\)"
        with pytest.raises(ValueError, match=match):
            scc.project_SCC(filter_by=[True, True])  # wrong number of lipids
            
        with pytest.raises(ValueError, match=match):
            scc.project_SCC(filter_by=[[True, True]])  # wrong number of frames

    def test_bins(self, scc):
        
        bins = np.linspace(0, 387, 388)
        scc_projection = scc.project_SCC(bins=bins)
        
        reference = {
            'extent': (-0.5, 386.5, -0.5, 386.5)
        }
        
        assert scc_projection.ax.images[0].get_extent() == reference['extent']
