
import pytest
import numpy as np
import MDAnalysis

from numpy.testing._private.utils import assert_array_almost_equal

from lipyphilic._simple_systems.simple_systems import (
    ONE_CHOL, ONE_CHOL_TRAJ)
from lipyphilic.lib.lateral_diffusion import MSD
 
 
class TestMSD:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    kwargs = {
        'lipid_sel': 'all'
    }
    
    def test_msd(self, universe):
        
        msd = MSD(universe, **self.kwargs)
        msd.run(stop=2)
    
        reference = {
            'n_residues': 1,
            'n_frames': 2,
            'msd': [[0.0, 0.04394855]],
            'lagtimes': [0.0, 5000.0]
        }
        
        assert msd.msd.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(msd.msd, reference['msd'])
        assert_array_almost_equal(msd.lagtimes, reference['lagtimes'])

    def test_msd_com_removal(self, universe):
        
        msd = MSD(universe, **self.kwargs, com_removal_sel="all")
        msd.run(stop=2)
        
        reference = {
            'n_residues': 1,
            'n_frames': 2,
            'msd': [[0.0, 0.0]],
            'lagtimes': [0.0, 5000.0]
        }
        
        assert msd.msd.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(msd.msd, reference['msd'])
        assert_array_almost_equal(msd.lagtimes, reference['lagtimes'])

    def test_diffusion_coefficient(self, universe):
        
        msd = MSD(universe, **self.kwargs, com_removal_sel="all")
        
        msd.msd = np.asarray([np.arange(100)])
        msd.lagtimes = np.arange(100)
        
        reference = {
            'd': 0.25 * 1e-5
        }
        
        d, sem = msd.diffusion_coefficient()
        
        assert d == reference['d']
        assert np.isnan(sem)  # we only have 1 molecule, so there is no SEM

    def test_diffusion_coefficient_start_stop(self, universe):
        
        msd = MSD(universe, **self.kwargs)
        
        msd.msd = np.asarray([np.arange(100)])
        msd.lagtimes = np.arange(100)
        
        reference = {
            'd': 0.25 * 1e-5
        }
        
        d, sem = msd.diffusion_coefficient(start_fit=20, stop_fit=80)
        
        assert d == reference['d']
        assert np.isnan(sem)  # we only have 1 molecule, so there is no SEM

    def test_diffusion_coefficient_lipid_sel(self, universe):
        
        msd = MSD(universe, **self.kwargs)
        
        msd.msd = np.asarray([np.arange(100)])
        msd.lagtimes = np.arange(100)
        
        reference = {
            'd': 0.25 * 1e-5
        }
        
        d, sem = msd.diffusion_coefficient(lipid_sel="all")
        
        assert d == reference['d']
        assert np.isnan(sem)  # we only have 1 molecule, so there is no SEM
