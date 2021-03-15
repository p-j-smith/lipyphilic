
import pytest
import numpy as np
import matplotlib
import MDAnalysis

from numpy.testing import assert_array_almost_equal
from lipyphilic.lib.plotting import ProjectionPlot

from lipyphilic._simple_systems.simple_systems import (
    HEX_LAT, HEX_LAT_BUMP_MID_MOL, HEX_LAT_OVERLAP)
from lipyphilic.lib.area_per_lipid import AreaPerLipid

matplotlib.use("Agg")
 
 
class TestAreaPerLipid:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    # first 50 residues are in the upper leaflet (1)
    # final 50 residues are in the upper leaflet (-1)
    kwargs = {
        'lipid_sel': 'name L C',
        'leaflets': np.array([[1]] * 50 + [[-1]] * 50)
    }
    
    @pytest.fixture(scope='class')
    def areas(self, universe):
        areas = AreaPerLipid(universe, **self.kwargs)
        areas.run()
        return areas

    def test_area_per_lipid(self, areas):
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'area': 200  # all lipids have the same area per lipid as they are on a hexagonal lattice
        }
    
        assert areas.areas.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(areas.areas, 200.0, decimal=8)
        
        
class TestAreaPerLipidOverlapping:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_OVERLAP)

    # first 50 residues are in the upper leaflet (1)
    # final 50 residues are in the upper leaflet (-1)
    kwargs = {
        'lipid_sel': 'name L C',
        'leaflets': np.array([[1]] * 50 + [[-1]] * 50)
    }
    
    @pytest.fixture(scope='class')
    def areas(self, universe):
        areas = AreaPerLipid(universe, **self.kwargs)
        areas.run()
        return areas
    
    def test_area_per_lipid(self, areas):
        
        # only the areas of 6 residues should be affected by these two overlapping atoms
        assert np.isclose(areas.areas, 200).sum() == 94


class TestAreaPerLipidMidplaneMol:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP_MID_MOL)

    # first 50 residues are in the upper leaflet (1)
    # final 50 residues are in the upper leaflet (-1), except:
    # CHOL78 is in the midplane
    leaflets = np.array([[1]] * 50 + [[-1]] * 50)
    leaflets[78] = 0
    kwargs = {
        'lipid_sel': 'name L C',
        'leaflets': leaflets
    }
    
    @pytest.fixture(scope='class')
    def areas(self, universe):
        areas = AreaPerLipid(universe, **self.kwargs)
        areas.run()
        return areas

    def test_area_per_lipid(self, areas):
    
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'min_area': 200  # all lipids have at least this area
        }
    
        assert areas.areas.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(np.nanmin(areas.areas), 200.0, decimal=8)
        assert np.isnan(areas.areas[78])
        assert sum(np.isnan(areas.areas)) == 1


class TestAreaPerLipidExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    def test_Exceptions(self, universe):
            
        match = "'leaflets' must either be a 1D array containing non-changing "
        with pytest.raises(ValueError, match=match):
            AreaPerLipid(
                universe=universe,
                lipid_sel="name L C",
                leaflets=None
            )
        
        with pytest.raises(ValueError, match=match):
            AreaPerLipid(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[[], []], [[], []]])  # cannot pass a 3D array
            )
        
        match = ("The shape of 'leaflets' must be \\(n_residues,\\), but 'lipid_sel' "
                 "generates an AtomGroup with 100 residues"
                 " and 'leaflets' has shape \\(99, 1\\).")
        with pytest.raises(ValueError, match=match):
            AreaPerLipid(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[1]] * 50 + [[-1]] * 49)  # one residue too few
            )
            
        match = ("The frames to analyse must be identical to those used "
                 "in assigning lipids to leaflets.")
        with pytest.raises(ValueError, match=match):
            areas = AreaPerLipid(
                universe=universe,
                lipid_sel="name L C",
                leaflets=np.array([[1, 1]] * 50 + [[-1, -1]] * 50)  # leaflets has two frames, apl one
            )
            areas.run()
            
            
class TestProjectArea:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)
    
    # first 50 residues are in the upper leaflet (1)
    # final 50 residues are in the upper leaflet (-1)
    kwargs = {
        'lipid_sel': 'name L C',
        'leaflets': np.array([[1]] * 50 + [[-1]] * 50)
    }
    
    @pytest.fixture(scope='class')
    def areas(self, universe):
        areas = AreaPerLipid(universe, **self.kwargs)
        areas.run()
        return areas
    
    def test_project_areas(self, areas):
        
        area_projection = areas.project_area()
        
        assert isinstance(area_projection, ProjectionPlot)
        assert_array_almost_equal(area_projection.values, areas.areas.mean())

    def test_filter_by(self, areas):
        
        area_projection = areas.project_area(filter_by=np.full(100, fill_value=True))
        
        assert_array_almost_equal(area_projection.values, areas.areas.mean())
        
    def test_filter_by_2D(self, areas):
        
        area_projection = areas.project_area(filter_by=np.full((100, 1), fill_value=True))
        
        assert_array_almost_equal(area_projection.values, areas.areas.mean())
        
    def test_filter_by_exception(self, areas):
        
        match = "The shape of `filter_by` must either be \(n_lipids, n_frames\) or \(n_lipids\)"
        with pytest.raises(ValueError, match=match):
            areas.project_area(filter_by=[True, True])  # wrong number of lipids
            
        with pytest.raises(ValueError, match=match):
            areas.project_area(filter_by=[[True, True]])  # wrong number of frames

    def test_bins(self, areas):
        
        bins = np.linspace(0, 100, 101)
        area_projection = areas.project_area(bins=bins)
        
        reference = {
            'extent': (-0.5, 99.5, -0.5, 99.5)
        }
        
        assert area_projection.ax.images[0].get_extent() == reference['extent']
