
import pytest
import numpy as np
import MDAnalysis

from numpy.testing._private.utils import assert_array_almost_equal, assert_array_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT_MONO
from lipyphilic.lib.z_thickness import ZThickness
 
 
class TestZThickness:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_MONO)

    kwargs = {
        'lipid_sel': 'name L C',
    }
    
    def test_ZThickness(self, universe):
        
        z_thickness = ZThickness(universe, **self.kwargs)
        z_thickness.run()
    
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'z_thickness': np.full((50, 1), fill_value=20)  # all lipids have a thickness of 20 Angstrom
        }
        
        assert z_thickness.z_thickness.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_equal(z_thickness.z_thickness, reference['z_thickness'])
        

class TestZThicknessAverage:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_MONO)
    
    @pytest.fixture(scope='class')
    def sn1_thickness(self, universe):
        sn1_thickness = ZThickness(universe, "name L")
        sn1_thickness.run()
        return sn1_thickness
    
    @pytest.fixture(scope='class')
    def sn2_thickness(self, universe):
        sn2_thickness = ZThickness(universe, "name C")
        sn2_thickness.run()
        return sn2_thickness
    
    def test_ZThickness_average(self, sn1_thickness):
        
        thickness = ZThickness.average(sn1_thickness, sn1_thickness)
        
        reference = {
            'n_residues': 25,
            'n_frames': 1,
            'z_thickness': np.full((25, 1), fill_value=20)
        }

        assert isinstance(thickness, ZThickness)
        assert thickness.z_thickness.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(thickness.z_thickness, reference['z_thickness'])
        
    def test_ZThickness_average_different_tails(self, sn1_thickness, sn2_thickness):
        
        thickness = ZThickness.average(sn1_thickness, sn2_thickness)
        
        reference = {
            'n_residues': 50,
            'n_frames': 1,
            'z_thickness': np.full((50, 1), fill_value=20)
        }

        assert thickness.z_thickness.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(thickness.z_thickness, reference['z_thickness'])
        

class TestZThicknessAverageExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_MONO)

    kwargs = {
        'lipid_sel': 'name L C',
    }
    
    @pytest.fixture(scope='class')
    def sn1_thickness(self, universe):
        sn1_thickness = ZThickness(universe, **self.kwargs)
        sn1_thickness.run()
        return sn1_thickness
    
    def test_Exceptions(self, universe, sn1_thickness):
            
        match = "sn1_thickness and sn2_thickness must have been run with the same frames"
        with pytest.raises(ValueError, match=match):
            sn2_thickness = ZThickness(
                universe=universe,
                **self.kwargs,
            )
            sn2_thickness.run(stop=0)
            ZThickness.average(sn1_thickness, sn2_thickness)
            
        with pytest.raises(ValueError, match=match):
            sn2_thickness = ZThickness(
                universe=universe,
                **self.kwargs,
            )
            sn2_thickness.run(stop=1)
            sn2_thickness.frames = np.array([10])
            ZThickness.average(sn1_thickness, sn2_thickness)
