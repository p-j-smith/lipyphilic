
import pytest
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from numpy.testing._private.utils import assert_array_almost_equal  # noqa: E402

from lipyphilic.lib.plotting import JointDensity  # noqa: E402
from lipyphilic.lib.plotting import ProjectionPlot  # noqa: E402
 
 
class TestProjectionPlot:
    
    kwargs = {
        'x_pos': [378.13067169],
        'y_pos': [62.28950091],
        'values': [0.12715311],
        'bins': np.linspace(0, 387, 388)
    }
    
    @pytest.fixture(scope="class")
    def projection(self):
        projection = ProjectionPlot(
            self.kwargs['x_pos'],
            self.kwargs['y_pos'],
            self.kwargs['values']
        )
        return projection
    
    def test_project_values(self, projection):
        
        projection.project_values(bins=self.kwargs['bins'])
        
        reference = {
            'n_values': 1,  # only one non-NaN value
            'values': 0.12715311
        }
        
        assert np.sum(~np.isnan(projection.statistic)) == reference['n_values']
        assert projection.statistic[~np.isnan(projection.statistic)] == reference['values']
        
    def test_interpolate(self, projection):
        
        projection.statistic = np.array(
            [
                [np.NaN, 0, np.NaN],
                [1, np.NaN, 1],
                [np.NaN, 0, np.NaN]
            ]
        )
        
        projection.interpolate(method="linear")
        
        reference = {
            'statistic': np.array(
                [
                    [1, 0.0, 0.5],
                    [1, 1, 1],
                    [1, 0.0, 0.5]
                ]
            )
        }
        
        assert_array_almost_equal(projection.statistic, reference['statistic'])
        
    def test_interpolate_no_tile(self, projection):
        
        projection.statistic = np.array(
            [
                [np.NaN, 0, np.NaN],
                [1, np.NaN, 1],
                [np.NaN, 0, np.NaN]
            ]
        )
        
        projection.interpolate(method="linear", tile=False)
        
        reference = {
            'statistic': np.array(
                [
                    [np.NaN, 0, np.NaN],
                    [1, 1, 1],
                    [np.NaN, 0, np.NaN]
                ]
            )
        }
        
        assert_array_almost_equal(projection.statistic, reference['statistic'])
    
    @pytest.fixture(scope="class")
    def projection_data(self):
        """A ProjectionPlot instance with the values calculated and interpolated.
        """
        
        projection = ProjectionPlot(
            self.kwargs['x_pos'],
            self.kwargs['y_pos'],
            self.kwargs['values']
        )
        projection.project_values(bins=self.kwargs['bins'])
        projection.interpolate(method="linear")
        
        return projection
    
    def test_plot_projection(self, projection_data):
        
        projection_data.plot_projection()
        
        reference = {
            'extent': (-0.5, 386.5, -0.5, 386.5)
        }
        
        assert projection_data.ax.images[0].get_extent() == reference['extent']
        assert isinstance(projection_data.cbar, matplotlib.colorbar.Colorbar)
    
    def test_plot_projection_existing_axis(self, projection_data):
        
        _, ax = plt.subplots(1)
        
        projection_data.plot_projection(ax=ax)
        
        assert projection_data.ax is ax
        
    def test_plot_projection_vmin_vmax(self, projection_data):
        
        projection_data.plot_projection(vmin=0.0, vmax=1.0)
        
        reference = {
            'cbar-ticks': np.linspace(0, 1, 6)
        }
        
        assert_array_almost_equal(projection_data.cbar.get_ticks(), reference['cbar-ticks'])
        
    def test_cmap(self, projection_data):
        
        projection_data.plot_projection(
            cmap="RdBu"
        )
        
        reference = {
            'cmap-name': 'RdBu'
        }
        
        assert projection_data._imshow.get_cmap().name == reference['cmap-name']
        
    def test_no_cbar(self, projection):
        
        # `projection_data` has already been used to create a plot with a cbar, and because the
        # scope of this object is this class, the cbar attribute is not None
        # So we use create a fresh plot with `projection`, which has not yet been used for plotting
        projection.plot_projection(cbar=False)
        
        assert projection.cbar is None


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


class TestInterpolate:
    
    @pytest.fixture(scope="class")
    def density(self):
        
        density = JointDensity(
            ob1=[],
            ob2=[]
        )
        
        # mock data without calling density.calc_density_2D()
        density.joint_mesh_values = np.array(
            [
                [np.NaN, 0, np.NaN],
                [1, np.NaN, 1],
                [np.NaN, 0, np.NaN]
            ]
        )
        
        return density
    
    def test_interpolate(self, density):
        
        density.interpolate(method="linear")
        
        reference = {
            'density': np.array(
                [
                    [0, 0, 0],
                    [1, 1, 1],
                    [0, 0, 0]
                ]
            )
        }
        
        assert np.sum(np.isnan(density.joint_mesh_values)) == 0
        assert_array_almost_equal(density.joint_mesh_values, reference['density'])
    
    def test_interpolate_fill_value(self, density):
        
        # Becayse density is a fixture with scope=class,
        # we need to reset the NaN values before interpolating
        density.joint_mesh_values = np.array(
            [
                [np.NaN, 0, np.NaN],
                [1, np.NaN, 1],
                [np.NaN, 0, np.NaN]
            ]
        )
        
        density.interpolate(method="linear", fill_value=100)
        
        reference = {
            'density': np.array(
                [
                    [100, 0, 100],
                    [1, 1, 1],
                    [100, 0, 100]
                ]
            )
        }
        
        assert np.sum(np.isnan(density.joint_mesh_values)) == 0
        assert_array_almost_equal(density.joint_mesh_values, reference['density'])
        
    def test_interpolate_PMF(self, density):
        
        # Becayse density is a fixture with scope=class,
        # we need to reset the NaN values before interpolating
        density.joint_mesh_values = np.array(
            [
                [np.NaN, 0, np.NaN],
                [1, np.NaN, 1],
                [np.NaN, 0, np.NaN]
            ]
        )
        
        density.temperature = 300  # The fill values should now be set to np.nanmax(density.joint_mesh_values), i.e 1.0
        density.interpolate(method="linear")
        
        reference = {
            'density': np.array(
                [
                    [1, 0, 1],
                    [1, 1, 1],
                    [1, 0, 1]
                ]
            )
        }
        
        assert np.sum(np.isnan(density.joint_mesh_values)) == 0
        assert_array_almost_equal(density.joint_mesh_values, reference['density'])
    

class TestPlotDensity:
    
    @pytest.fixture(scope="class")
    def density(self):
        
        density = JointDensity(
            ob1=[],
            ob2=[]
        )
        
        # mock data without calling density.calc_density_2D()
        # Below is the data from running calc_density_2D on
        # the ob1 and ob2 data in TestJointDensity
        density.joint_mesh_values = np.zeros((180, 20), dtype=float)
        
        non_zero_angles = np.array(
            [
                5, 7, 14, 21, 26, 27, 27, 31, 47, 87,
                88, 115, 140, 153, 155, 157, 159, 162, 164, 165,
                166, 171, 173, 176
            ]
        )
        non_zero_heights = np.array(
            [
                16, 16, 16, 16, 15, 12, 16, 14, 12, 10,
                10, 8, 5, 3, 4, 3, 3, 5, 3, 3,
                3, 3, 3, 3
            ]
        )
        non_zero_values = np.array(
            [
                0.04, 0.04, 0.04, 0.04, 0.04,
                0.04, 0.04, 0.04, 0.04, 0.04,
                0.04, 0.04, 0.04, 0.04, 0.04,
                0.04, 0.04, 0.04, 0.08, 0.04,
                0.04, 0.04, 0.04, 0.04
            ]
        )
        
        density.joint_mesh_values[non_zero_angles, non_zero_heights] = non_zero_values
        
        density.ob1_mesh_bins, density.ob2_mesh_bins = np.meshgrid(
            np.linspace(0.5, 179.5, 180),
            np.linspace(-9.5, 9.5, 20)
        )
        
        return density

    def test_plot_density(self, density):
        
        density.plot_density()
        
        reference = {
            'contour-vertices': np.array(
                [
                    [153.5, -7.0],
                    [154.0, -6.5],
                    [153.5, -6.0],
                    [153.0, -6.5],
                    [153.5, -7.0]
                ]
            ),
            'cbar-orientation': 'vertical',
            'cbar-ylabel': 'Probability density',
            'cbar-aspect': 30,
            'cbar-x-extent': (0.78375, 0.809417),
            'cbar-ticks': np.linspace(0, 0.08, 9)
        }
        
        assert_array_almost_equal(density.ax.collections[1].get_paths()[0].vertices.round(1), reference['contour-vertices'])
        assert density.cbar.orientation == reference['cbar-orientation']
        assert density.cbar.ax.get_ylabel() == reference['cbar-ylabel']
        assert_array_almost_equal((density.cbar.ax.get_position().x0, density.cbar.ax.get_position().x1), reference['cbar-x-extent'], decimal=1)
        assert_array_almost_equal(density.cbar.get_ticks(), reference['cbar-ticks'])

    def test_plot_density_existing_axis(self, density):
        
        _, ax = plt.subplots(1)
        
        density.plot_density(ax=ax)
        
        assert density.ax is ax
        
    def test_vmin_vmax(self, density):
        
        density.plot_density(vmin=0.0, vmax=1.0)
        
        reference = {
            'cbar-ticks': np.linspace(0, 1, 6)
        }
        
        assert_array_almost_equal(density.cbar.get_ticks(), reference['cbar-ticks'])
            
    def test_labels(self, density):
        
        density.plot_density(
            title="My Title",
            xlabel="My x-label",
            ylabel="My y-label"
        )
        
        reference = {
            'title': 'My Title',
            'xlabel': 'My x-label',
            'ylabel': 'My y-label'
        }
        
        assert density.ax.get_title(loc="left") == reference['title']
        assert density.ax.get_xlabel() == reference['xlabel']
        assert density.ax.get_ylabel() == reference['ylabel']
        
    def test_cmap(self, density):
        
        density.plot_density(
            cmap="RdBu"
        )
        
        reference = {
            'cmap-name': 'RdBu'
        }
        
        assert density._imshow.get_cmap().name == reference['cmap-name']
        
    def test_cbar_kws(self, density):
        
        density.plot_density(
            cbar_kws={
                'label': 'My label',
                'aspect': 10,
                'pad': 0.1
            }
        )
        
        reference = {
            'label': 'My label',
            'aspect': 10,
            'x-extent': (0.78375, 0.86075),
        }
        
        assert density.cbar.ax.get_ylabel() == reference['label']
        assert density.cbar.ax.get_aspect() == reference['aspect']
        assert_array_almost_equal((density.cbar.ax.get_position().x0, density.cbar.ax.get_position().x1), reference['x-extent'], decimal=1)
        
    def test_contour_kws(self, density):
        
        density.plot_density(
            contour_kws={"colors": "red"}
        )
        
        reference = {
            'colors': 'red'
        }
        
        assert density._contours.colors == reference['colors']
        
    def test_clabel_kws(self, density):
        
        density.plot_density(
            contour_labels=[0, 1, 2, 3, 4],
            clabel_kws={"fontsize": 1}  # small font size requried otherwise no labels added to this plot
        )
        
        reference = {
            'n_labels': 43  # total number of labels added to the contour lines
        }
        
        assert len(density._clabels) == reference['n_labels']


class TestPlotPMF:
    
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
        density.calc_density_2D(
            bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']),
            temperature=300
            )
        density.interpolate(method="cubic")
        
        return density
    
    @pytest.fixture(scope="class")
    def density_310K(self):
        
        density = JointDensity(
            self.kwargs['ob1'],
            self.kwargs['ob2']
        )
        density.calc_density_2D(
            bins=(self.kwargs['angle-bins'], self.kwargs['height-bins']),
            temperature=310
            )
        density.interpolate()
        
        return density

    def test_plot_PMF(self, density):
        
        density.plot_density()
        
        reference = {
            'contour-vertices': np.array(
                [
                    [153.5, -7.0],
                    [154.0, -6.5],
                    [153.5, -6.0],
                    [153.0, -6.5],
                    [153.5, -7.0]
                ]
            ),
            'cbar-orientation': 'vertical',
            'cbar-ylabel': 'PMF',
            'cbar-aspect': 30,
            'cbar-x-extent': (0.78375, 0.809417),
            'cbar-ticks': np.linspace(0, 3, 7)
        }
        
        # assert_array_almost_equal(density.ax.collections[1].get_paths()[0].vertices.round(1), reference['contour-vertices'])
        assert density.cbar.orientation == reference['cbar-orientation']
        assert_array_almost_equal((density.cbar.ax.get_position().x0, density.cbar.ax.get_position().x1), reference['cbar-x-extent'], decimal=1)
        assert_array_almost_equal(density.cbar.get_ticks(), reference['cbar-ticks'])
    
    def test_PMF_difference(self, density, density_310K):
        
        density.plot_density(difference=density_310K)
        
        reference = {
            'cbar-ticks': np.linspace(-1, 1, 9)
        }
                
        assert_array_almost_equal(density.cbar.get_ticks(), reference['cbar-ticks'])
        
    def test_no_cbar(self, density_310K):
            
        density_310K.plot_density(cbar=False)
        
        assert density_310K.cbar is None
