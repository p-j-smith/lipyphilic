
import pytest
import numpy as np
import MDAnalysis

from numpy.testing import assert_array_equal

from lipyphilic._simple_systems.simple_systems import (
    HEX_LAT, HEX_LAT_BUMP)
from lipyphilic.lib.z_positions import ZPositions
 
 
class TestZPositions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C',
        'height_sel': 'name C',
        'n_bins': 1
    }
    
    @pytest.fixture(scope='class')
    def z_positions(self, universe):
        z_positions = ZPositions(universe, **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
    
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'z_positions': np.full((50, 1), fill_value=10)
        }
        
        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference['z_positions'][25:] = -10
        
        assert z_positions.z_positions.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_equal(z_positions.z_positions, reference['z_positions'])
        
        
class TestZPositionsOneAtom:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C',
        'n_bins': 1
    }
    
    @pytest.fixture(scope='class')
    def z_positions(self, universe):
        height_indices = " ".join(universe.atoms[2::4].indices.astype(str))  # one atom per cholesterol
        z_positions = ZPositions(universe, height_sel=f"index {height_indices}", **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
    
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'z_positions': np.full((50, 1), fill_value=10)
        }
        
        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference['z_positions'][25:] = -10
        
        assert z_positions.z_positions.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_equal(z_positions.z_positions, reference['z_positions'])


class TestZPositionsUndulating:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_BUMP)

    kwargs = {
        'lipid_sel': 'name L C',
        'height_sel': 'name C',
        'n_bins': 10  # with fewer bins the z positions will no longer be +/- 10 Angstrom
    }
    
    @pytest.fixture(scope='class')
    def z_positions(self, universe):
        z_positions = ZPositions(universe, **self.kwargs)
        z_positions.run()
        return z_positions

    def test_z_positions(self, z_positions):
    
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'z_positions': np.full((50, 1), fill_value=10)
        }
        
        # the lower leaflet (final 25 residues) cholesterols are at -10 Angstrom
        reference['z_positions'][25:] = -10
        
        assert z_positions.z_positions.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_equal(z_positions.z_positions, reference['z_positions'])
