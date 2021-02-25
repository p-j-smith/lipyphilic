
import pytest
import numpy as np
import MDAnalysis

from numpy.testing import assert_array_almost_equal

from lipyphilic._simple_systems.simple_systems import HEX_LAT
from lipyphilic.lib.assign_leaflets import AssignLeaflets
from lipyphilic.lib.registration import Registration
 
 
class TestRegistration:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)
    
    @pytest.fixture(scope='class')
    def leaflets(self, universe):
        leaflets = AssignLeaflets(universe, lipid_sel="name L C")
        leaflets.run()
        return leaflets

    def test_registred(self, leaflets):
    
        registration = Registration(
            leaflets=leaflets,
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
        
    def test_antiregistred(self, leaflets):
        
        registration = Registration(
            leaflets=leaflets,
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

    def test_nbins100(self, leaflets):
        
        registration = Registration(
            leaflets=leaflets,
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
        
    def test_nbins1000_sd0(self, leaflets):
        
        # The bins are 0.1 Anstrom wide, and gaussian_sd=0 means there is no spead
        # in density
        # So, registration should be
        registration = Registration(
            leaflets=leaflets,
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
        
    def test_filter_by_registered(self, leaflets):
        
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        
        registration = Registration(
            leaflets=leaflets,
            filter_by=filter_by
        )
        registration.run()
    
        reference = {
            'n_frames': 1,
            'registration': 1
        }
    
        assert registration.registration.size == reference['n_frames']
        assert_array_almost_equal(registration.registration, reference['registration'], decimal=4)
        
    def test_filter_by_2D_registered(self, leaflets):
        
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        
        registration = Registration(
            leaflets=leaflets,
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
    
    @pytest.fixture(scope='class')
    def leaflets(self, universe):
        leaflets = AssignLeaflets(universe, lipid_sel="name L C")
        leaflets.run()
        return leaflets
    
    def test_Exceptions(self, leaflets):
            
        match = "leaflets must be of type AssignLeaflets"
        with pytest.raises(ValueError, match=match):
            Registration(
                leaflets=leaflets.leaflets  # this is a NumPy array
            )
            
        filter_by = np.zeros(100)
        filter_by[1::2] = 1
        match = "'filter_by' must either be a 1D array containing non-changing boolean"
        with pytest.raises(ValueError, match=match):
            Registration(
                leaflets=leaflets,
                filter_by=np.array(None)
            )
            
        match = "The shape of 'filter_by' must be \(n_residues,\)"  # noqa:W605
        with pytest.raises(ValueError, match=match):
            Registration(
                leaflets=leaflets,
                filter_by=filter_by[:99]
            )
