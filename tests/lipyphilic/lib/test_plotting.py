
import pytest
import numpy as np
import scipy
import matplotlib

from numpy.testing._private.utils import assert_array_almost_equal

from lipyphilic.lib.plotting import JointDensity
 
matplotlib.use("Agg")
 
 
class TestJointDensity:
    
    # ZAngle of ONE_CHOL_TRAJ
    angle = np.array(
        [
            [
                164.07437819, 165.14768719, 155.39577416, 176.70383898, 159.41782711,
                153.8045642, 115.96163667, 26.42859516, 47.08837727, 140.84573092,
                88.29349009, 21.31531296, 27.73287278, 5.1672467, 27.89790106,
                31.60123903, 7.97378798, 14.0942026, 87.08317218, 164.64607065,
                171.55622297, 166.27983484, 162.90936118, 173.60002084, 157.20084899
            ]
        ]
    )
    # ZPosition of ONE_CHOL_TRAJ
    height = np.array(
        [
            [
                -6.13125229, -6.58125257, -5.81000042, -6.78375101, -6.68249846,
                -6.56250143, -1.96625185, 5.66625118, 2.03750038, -4.92875004,
                0.0, 6.65625191, 6.26000071, 6.77999878, 2.11999893,
                4.39666557, 6.79874945, 6.61624622, 0.67749882, -6.49000096,
                -6.88750029, -6.85375118, -4.04749966, -6.94000149, -6.5612483
            ]
        ]
    )
    
    kwargs = {
        'ob1': angle,
        'ob2': height,
        'angle-bins': np.linspace(0, 180, 181),
        'height-bins': np.linspace(-10, 10, 21),
        'temperature': 300
    }
    
    @pytest.fixture(scope="class")
    def density(self):
        density = JointDensity(
            self.kwargs['ob1'],
            self.kwargs['ob2']
        )
        return density
    
    def test_joint_density_exceptions(self):
            
        match = "`ob1` and `ob2` must be arrays of the same shape."
        with pytest.raises(ValueError, match=match):
            JointDensity(
                self.kwargs['ob1'],
                self.kwargs['ob2'][:-1]
            )
    
    def test_calc_density(self, density):
        
        density.calc_density_2D(bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']))
        
        reference = {
            'range': [0.0, 0.04, 0.08]  # min density, min non-zero density, max density
        }
        
        actual = {
            'range': [
                np.nanmin(density.joint_mesh_values),
                np.min(density.joint_mesh_values[density.joint_mesh_values > 0.0]),
                np.nanmax(density.joint_mesh_values)
            ]
        }
        
        assert_array_almost_equal(actual['range'], reference['range'])
        
    def test_calc_density_filter(self, density):
        
        density.calc_density_2D(
            bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']),
            filter_by=self.kwargs['ob2'] > 0.0
        )
        
        reference = {
            'range': [0.0, 0.1, 0.1]  # min density, min non-zero density, max density
        }
        
        actual = {
            'range': [
                np.nanmin(density.joint_mesh_values),
                np.min(density.joint_mesh_values[density.joint_mesh_values > 0.0]),
                np.nanmax(density.joint_mesh_values)
            ]
        }
        
        assert_array_almost_equal(actual['range'], reference['range'])
        
    def test_calc_density_exceptions(self, density):
        
        match = "`filter_by` must be an array with the same shape as `ob1` and `ob2`."
        with pytest.raises(ValueError, match=match):
            density.calc_density_2D(
                bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']),
                filter_by=self.kwargs['ob2'][:, :-1] > 0.0  # mask is too small in dimension 1
            )
        
    def test_calc_PMF(self, density):
        
        density.calc_density_2D(
            bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']),
            temperature=self.kwargs['temperature']
            )
        
        density_range = [0.08, 0.04]  # max density, min non-zero density
        pmf_range = -(scipy.constants.Boltzmann * scipy.constants.Avogadro / 4184) * self.kwargs['temperature'] * np.log(density_range)
        
        reference = {
            'range': pmf_range
        }
        
        actual = {
            'range': [
                np.min(density.joint_mesh_values[density.joint_mesh_values > 0.0]),
                np.nanmax(density.joint_mesh_values)
            ]
        }
        
        assert_array_almost_equal(actual['range'], reference['range'])
            