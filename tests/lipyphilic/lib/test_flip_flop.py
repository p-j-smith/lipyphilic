
import pytest
import numpy as np
import MDAnalysis

from numpy.testing._private.utils import assert_array_equal

from lipyphilic._simple_systems.simple_systems import (
    ONE_CHOL, ONE_CHOL_TRAJ)
from lipyphilic.lib.flip_flop import FlipFlop
 
  
class TestFlipFlop:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    @staticmethod
    @pytest.fixture(scope='class')
    def leaflets():
        """Leaflets the ONE_CHOL molecule belongs to at each frame.
        """
        leaflets = np.array(
            [
                -1, -1, -1, -1, -1, -1, 0, 0, 1, 0,
                0, 1, 1, 1, 1, 1, 1, 1, 0, -1,
                -1, -1, -1, -1, -1
            ],
            dtype=np.int8
        )
        #  We need to make sure leaflets has the correct shape
        #  The shape is (n_residues=1, n_frames=25)
        leaflets = leaflets[np.newaxis, :]
        return leaflets

    def test_flip_flop(self, universe, leaflets, lipid_sel="name ROH"):
        
        # 3 flip-flop events with the second one a failure (remains in the upper leaflet)
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets,
            frame_cutoff=1,
        )
        flip_flop.run()
        
        reference = {
            'events': [
                [0, 5, 8, 1],
                [0, 8, 11, 1],
                [0, 17, 19, -1]
            ],
            'success': ["Success", "Fail", "Success"]
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])
        
    def test_flip_flop_framecut2(self, universe, leaflets, lipid_sel="name ROH"):
            
        # only two flip-flop events found now, both successful
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets,
            frame_cutoff=2,
        )
        flip_flop.run()
        
        reference = {
            'events': [
                [0, 5, 11, 1],
                [0, 17, 19, -1]
            ],
            'success': ["Success", "Success"]
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])

    def test_flip_flop_framecut8(self, universe, leaflets, lipid_sel="name ROH"):
                
        # No flip-flop events found now
        # Cholesterol doesn't remain in opposing lealet for long enough
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets,
            frame_cutoff=8,
        )
        flip_flop.run()
        
        reference = {
            'events': np.array([[], [], [], []]).T,
            'success': np.array([])
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])
        
    def test_flip_flop_same_leaflet(self, universe, leaflets, lipid_sel="name ROH"):
                
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=np.ones_like(leaflets, dtype=np.int8),
        )
        flip_flop.run()
        
        reference = {
            'events': np.array([[], [], [], []]).T,
            'success': np.array([])
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])
        
    def test_flip_flop_upper_mid_only(self, universe, leaflets, lipid_sel="name ROH"):
        
        leaflets = np.ones_like(leaflets, dtype=np.int8)
        leaflets[0, 1] = 0
        
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets
        )
        flip_flop.run()
        
        reference = {
            'events': np.array([[0], [0], [2], [1]]).T,
            'success': np.array(["Fail"])
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])
        
    def test_flip_flop_lower_mid_only(self, universe, leaflets, lipid_sel="name ROH"):
        
        leaflets = np.ones_like(leaflets, dtype=np.int8) * -1
        leaflets[0, 1] = 0
        
        flip_flop = FlipFlop(
            universe=universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets
        )
        flip_flop.run()
        
        reference = {
            'events': np.array([[0], [0], [2], [-1]]).T,
            'success': np.array(["Fail"])
        }
        
        assert_array_equal(flip_flop.flip_flops, reference['events'])
        assert_array_equal(flip_flop.flip_flop_success, reference['success'])


class TestAreaPerLipidExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(ONE_CHOL, ONE_CHOL_TRAJ)

    @staticmethod
    @pytest.fixture(scope='class')
    def leaflets():
        """Leaflets the ONE_CHOL molecule belongs to at each frame.
        """
        leaflets = np.array(
            [
                -1, -1, -1, -1, -1, -1, 0, 0, 1, 0,
                0, 1, 1, 1, 1, 1, 1, 1, 0, -1,
                -1, -1, -1, -1, -1
            ],
            dtype=np.int8
        )
        #  We need to make sure leaflets has he correct shape
        #  The shape is (n_residues=1, n_frames=25)
        leaflets = leaflets[np.newaxis, :]
        return leaflets

    def test_Exceptions(self, universe, leaflets, lipid_sel="name ROH"):
            
        match = ("'leaflets' must be a 2D array of shape \\(n_residues, n_frames\\)"
                 " containing the leaflet id of each lipid at each frame."
                 )
        with pytest.raises(ValueError, match=match):
            FlipFlop(
                universe=universe,
                lipid_sel=lipid_sel,
                leaflets=leaflets[0]  # wrong dimensions
            )
        with pytest.raises(ValueError, match=match):
            FlipFlop(
                universe=universe,
                lipid_sel=lipid_sel,
                leaflets=np.concatenate((leaflets, leaflets))  # wrong number of residues
            )
            
        match = ("The frames to analyse must be identical to those used "
                 "in assigning lipids to leaflets."
                 )
        with pytest.raises(ValueError, match=match):
            flip_flop = FlipFlop(
                universe=universe,
                lipid_sel=lipid_sel,
                leaflets=leaflets[:, ::2]  # not enough frames for the trajectory
            )
            flip_flop.run(
                start=None,
                stop=None,
                step=None
            )
        with pytest.raises(ValueError, match=match):
            flip_flop = FlipFlop(
                universe=universe,
                lipid_sel=lipid_sel,
                leaflets=leaflets
            )
            flip_flop.run(
                start=None,
                stop=None,
                step=2  # different number of frames to what is in 'leaflets'
            )
            
        match = "'frame_cutoff' must be greater than or equal to 1"
        with pytest.raises(ValueError, match=match):
            FlipFlop(
                universe=universe,
                lipid_sel=lipid_sel,
                leaflets=leaflets,
                frame_cutoff=0
            )
