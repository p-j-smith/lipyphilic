
import pathlib
import pytest
import numpy as np
import MDAnalysis

from numpy.testing import assert_array_almost_equal, assert_raises

from lipyphilic._simple_systems.simple_systems import (
    HEX_LAT_TRANS, HEX_LAT_TRANS_TRAJ, HEX_LAT_SPLIT_Z
)
from lipyphilic.transformations import nojump, center_membrane


class TestNoJump:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_TRANS, HEX_LAT_TRANS_TRAJ)

    # All atoms were translated in y by 5 Angstrom at each frame
    # Upper leaflet atoms move in the negative  direction
    # Lower leaflet move in the positive y direction
    distance_moved = 5
    n_frames = 3
    n_atoms = 200
    
    expected_y_diffs = np.full((n_frames - 1, n_atoms), fill_value=distance_moved)
    expected_y_diffs[:, :100] *= -1

    reference = {
        "x_diffs": np.zeros((n_frames - 1, n_atoms)),
        "y_diffs": expected_y_diffs,
        "z_diffs": np.zeros((n_frames - 1, n_atoms)),
    }
    
    def test_without_no_jump(self, universe):
        
        atoms = universe.select_atoms("all")
        
        # Get positions at each frame
        all_positions = np.zeros((atoms.n_atoms, universe.trajectory.n_frames, 3))
        for ts in universe.trajectory:
            all_positions[:, ts.frame] = atoms.positions
            
        # And the distance they move at each frame
        x_diffs, y_diffs, z_diffs = np.diff(all_positions, axis=1).T

        # Atoms are stationary in x, but get wrapped when the box shrinks
        assert_raises(
            AssertionError,
            assert_array_almost_equal,
            self.reference["x_diffs"], x_diffs, decimal=5
        )
        
        # Atoms are translated in y
        assert_raises(
            AssertionError,
            assert_array_almost_equal,
            self.reference["y_diffs"], y_diffs, decimal=5
        )
        
        # Atoms are not translated in z and the box does no shrink enough to cause wrapping in z
        np.testing.assert_array_almost_equal(self.reference["z_diffs"], z_diffs, decimal=5)
        
    def test_with_no_jump(self, universe):
        
        atoms = universe.select_atoms("all")
        universe.trajectory.add_transformations(nojump(atoms, nojump_z=True))
        
        # Get positions at each frame
        all_positions = np.zeros((atoms.n_atoms, universe.trajectory.n_frames, 3))
        for ts in universe.trajectory:
            all_positions[:, ts.frame] = atoms.positions
            
        # And the distance they move at each frame
        x_diffs, y_diffs, z_diffs = np.diff(all_positions, axis=1).T

        np.testing.assert_array_almost_equal(self.reference["x_diffs"], x_diffs, decimal=5)
        np.testing.assert_array_almost_equal(self.reference["y_diffs"], y_diffs, decimal=5)
        np.testing.assert_array_almost_equal(self.reference["z_diffs"], z_diffs, decimal=5)

 
class TestNoJumpStatic:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        
        # Write new nojump trajectory
        u = MDAnalysis.Universe(HEX_LAT_TRANS, HEX_LAT_TRANS_TRAJ)
        
        xtc_folder = pathlib.Path(HEX_LAT_TRANS_TRAJ).parent
        filename = "_HexGrid-2AtomsPerLipid-TranslatedIn_y_nojump.xtc"
        HEX_LAT_TRANS_TRAJ_NOJUMP = xtc_folder.joinpath(filename).as_posix()

        atoms = u.select_atoms("all")
        u.trajectory.add_transformations(
            nojump(
                atoms,
                nojump_z=True,
                filename=HEX_LAT_TRANS_TRAJ_NOJUMP
            )
        )
        
        return MDAnalysis.Universe(HEX_LAT_TRANS, HEX_LAT_TRANS_TRAJ_NOJUMP)
    
    # All atoms were translated in y by 5 Angstrom at each frame
    # Upper leaflet atoms move in the negative  direction
    # Lower leaflet move in the positive y direction
    distance_moved = 5
    n_frames = 3
    n_atoms = 200
    
    expected_y_diffs = np.full((n_frames - 1, n_atoms), fill_value=distance_moved)
    expected_y_diffs[:, :100] *= -1

    reference = {
        "x_diffs": np.zeros((n_frames - 1, n_atoms)),
        "y_diffs": expected_y_diffs,
        "z_diffs": np.zeros((n_frames - 1, n_atoms)),
    }
    
    def test_no_jump_static(self, universe):
        
        atoms = universe.select_atoms("all")
        
        # Get positions at each frame
        all_positions = np.zeros((atoms.n_atoms, universe.trajectory.n_frames, 3))
        for ts in universe.trajectory:
            all_positions[:, ts.frame] = atoms.positions
            
        # And the distance they move at each frame
        x_diffs, y_diffs, z_diffs = np.diff(all_positions, axis=1).T

        np.testing.assert_array_almost_equal(self.reference["x_diffs"], x_diffs, decimal=5)
        np.testing.assert_array_almost_equal(self.reference["y_diffs"], y_diffs, decimal=5)
        np.testing.assert_array_almost_equal(self.reference["z_diffs"], z_diffs, decimal=5)
    

class TestCenterMembrane:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT_SPLIT_Z)
    
    reference = {
        "upper_leaflet": np.arange(100),
        "lower_leaflet": np.arange(100, 200),
        "broken_bilayer_height": 92,
        "broken_bilayer_midpoint": 0,
        "correct_bilayer_height": 8,
        "bilayer_midpoint": 50
    }
    
    def test_no_center(self, universe):
        
        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])
        
        assert bilayer_height == self.reference["broken_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]
        
    def test_center_frame(self, universe):
        
        membrane = universe.select_atoms("name L C")
        ts = universe.trajectory[0]
        center_membrane(ag=membrane, shift=5)(ts)
        
        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])
        
        assert bilayer_height == self.reference["correct_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]
        
    def test_center_trajectory(self):
        
        universe = MDAnalysis.Universe(HEX_LAT_SPLIT_Z)
        
        membrane = universe.select_atoms("name L C")
        universe.trajectory.add_transformations(center_membrane(ag=membrane))
        
        upper_z_pos = universe.atoms.positions[self.reference["upper_leaflet"], 2]
        lower_z_pos = universe.atoms.positions[self.reference["lower_leaflet"], 2]
        bilayer_height = abs(np.mean(upper_z_pos) - np.mean(lower_z_pos))
        bilayer_midpoint = np.mean([upper_z_pos, lower_z_pos])
        
        assert bilayer_height == self.reference["correct_bilayer_height"]
        assert bilayer_midpoint == self.reference["bilayer_midpoint"]
