
import pytest
import numpy as np
import MDAnalysis

from numpy.testing import assert_array_almost_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT
from lipyphilic.lib.registration import Registration
 
 
class TestRegistration:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)
    
    kwargs = {
        'leaflets': np.array([[1]] * 50 + [[-1]] * 50)
    }

    def test_registred(self, universe):
    
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L C",
            lower_sel="name L C"
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'])
        
    def test_antiregistred(self, universe):
        
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L",
            lower_sel="name C"
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': -1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'])

    def test_nbins100(self, universe):
        
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L C",
            lower_sel="name L C",
            n_bins=100
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'])
        
    def test_nbins1000_sd0(self, universe):
        
        # The bins are 0.1 Anstrom wide, and gaussian_sd=0 means there is no spead
        # in density
        # So, registration should be almost 0
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L",
            lower_sel="name C",
            n_bins=1000,
            gaussian_sd=0
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 0
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'], decimal=4)
        
    def test_filter_by_registered(self, universe):
        
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L C",
            lower_sel="name L C",
            filter_by=filter_by
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'], decimal=4)
        
    def test_filter_by_2D_registered(self, universe):
        
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        
        registration = Registration(
            universe=universe,
            leaflets=self.kwargs['leaflets'],
            upper_sel="name L C",
            lower_sel="name L C",
            filter_by=filter_by[:, np.newaxis]
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'], decimal=4)


class TestRegistrationExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)
    
    kwargs = {
        'leaflets': np.array([[1]] * 50 + [[-1]] * 50)
    }
    
    def test_Exceptions(self, universe):
        
        match = "'leaflets' must either be a 1D array containing non-changing "
        with pytest.raises(ValueError, match=match):
            Registration(
                universe=universe,
                upper_sel="name L C",
                lower_sel="name L C",
                leaflets=None
            )
        
        with pytest.raises(ValueError, match=match):
            Registration(
                universe=universe,
                upper_sel="name L C",
                lower_sel="name L C",
                leaflets=np.array([[[], []], [[], []]])  # cannot pass a 3D array
            )
        
        match = ("The shape of 'leaflets' must be \\(n_residues,\\), but 'lipid_sel' "
                 "generates an AtomGroup with 100 residues"
                 " and 'leaflets' has shape \\(99, 1\\).")
        with pytest.raises(ValueError, match=match):
            Registration(
                universe=universe,
                upper_sel="name L C",
                lower_sel="name L C",
                leaflets=np.array([[1]] * 50 + [[-1]] * 49)  # one residue too few
            )
        
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        match = "'filter_by' must either be a 1D array containing non-changing boolean"
        with pytest.raises(ValueError, match=match):
            Registration(
                universe=universe,
                upper_sel="name L C",
                lower_sel="name L C",
                leaflets=self.kwargs['leaflets'],
                filter_by=np.array(None)
            )
            
        match = "The shape of 'filter_by' must be \(n_residues,\)"  # noqa:W605
        with pytest.raises(ValueError, match=match):
            Registration(
                universe=universe,
                upper_sel="name L C",
                lower_sel="name L C",
                leaflets=self.kwargs['leaflets'],
                filter_by=filter_by[:99]
            )
