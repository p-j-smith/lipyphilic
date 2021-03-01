
import pytest
import numpy as np
import MDAnalysis

from numpy.testing._private.utils import assert_array_almost_equal

from lipyphilic._simple_systems.simple_systems import (
    ONE_CHOL, ONE_CHOL_TRAJ)
from lipyphilic.lib.z_angles import ZAngles
 
 
class TestZAngles:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    kwargs = {
        'atom_A_sel': 'name R5',
        'atom_B_sel': 'name ROH',
    }
    
    def test_z_angles_degrees(self, universe):
        
        z_angles = ZAngles(universe, **self.kwargs, rad=False)
        z_angles.run()
    
        reference = {
            'n_residues': 1,
            'n_frames': 25,
            'z_angles': np.array(
                [
                    [
                        15.92562181, 14.85231281, 24.60422584, 3.29616102, 20.58217289,
                        26.1954358, 64.03836333, 153.57140484, 132.91162273, 39.15426908,
                        91.70650991, 158.68468704, 152.26712722, 174.8327533, 152.10209894,
                        148.39876097, 172.02621202, 165.9057974, 92.91682782, 15.35392935,
                        8.44377703, 13.72016516, 17.09063882, 6.39997916, 22.79915101
                    ]
                ]
            )
        }
        
        assert z_angles.z_angles.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(z_angles.z_angles, reference['z_angles'])
        
    def test_z_angles_radians(self, universe):
        
        z_angles = ZAngles(universe, **self.kwargs, rad=True)
        z_angles.run()
    
        reference = {
            'n_residues': 1,
            'n_frames': 25,
            'z_angles': np.array(
                [
                    [
                        15.92562181, 14.85231281, 24.60422584, 3.29616102, 20.58217289,
                        26.1954358, 64.03836333, 153.57140484, 132.91162273, 39.15426908,
                        91.70650991, 158.68468704, 152.26712722, 174.8327533, 152.10209894,
                        148.39876097, 172.02621202, 165.9057974, 92.91682782, 15.35392935,
                        8.44377703, 13.72016516, 17.09063882, 6.39997916, 22.79915101
                    ]
                ]
            )
        }
        
        assert z_angles.z_angles.shape == (reference['n_residues'], reference['n_frames'])
        assert_array_almost_equal(z_angles.z_angles, np.deg2rad(reference['z_angles']))


class TestZAnglesExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    kwargs = {
        'atom_A_sel': 'name R5',
        'atom_B_sel': 'name PO4',  # this bead doesn't exist in CHOL
    }

    def test_Exceptions(self, universe):
            
        match = "atom_A_sel and atom_B_sel must select the same number of atoms"
        with pytest.raises(ValueError, match=match):
            ZAngles(
                universe=universe,
                **self.kwargs
            )
