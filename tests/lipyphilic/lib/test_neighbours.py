

import pytest
import numpy as np
from numpy.testing import assert_array_equal
import MDAnalysis
from MDAnalysis.exceptions import NoDataError

from lipyphilic._simple_systems.simple_systems import HEX_LAT
from lipyphilic.lib.neighbours import Neighbours
 
 
class TestNeighbours:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C'
    }
    
    def test_neighbours_cutoff12(self, universe):
    
        neighbours = Neighbours(universe, **self.kwargs, cutoff=12)
        neighbours.run()
    
        # it's a hexagonal lattice
        # with cutoff=12, every lipid should have 6 neighbours
        reference = {
            'n_frames': 1,
            'n_residues': 100,  # there are 200 atoms but 100 lipids in total
            'n_neighbours': 6,
        }
    
        assert neighbours.neighbours.shape[0] == (reference['n_frames'])
        assert neighbours.neighbours[0].shape == (reference['n_residues'], reference['n_residues'])
        assert (np.sum(neighbours.neighbours[0].toarray(), axis=0) == reference['n_neighbours']).all()
        
    def test_neighbours_cutoff10(self, universe):
        
        neighbours = Neighbours(universe, **self.kwargs, cutoff=10)
        neighbours.run()
    
        # it's a hexagonal lattice, but each hexagon is irregular (two sides longer than the other 4)
        # with cutoff=10, every lipid should have 2 neighbours
        reference = {
            'n_frames': 1,
            'n_residues': 100,  # there are 200 atoms but 100 lipids in total
            'n_neighbours': 2,
        }
    
        assert neighbours.neighbours.shape[0] == (reference['n_frames'])
        assert neighbours.neighbours[0].shape == (reference['n_residues'], reference['n_residues'])
        assert (np.sum(neighbours.neighbours[0].toarray(), axis=0) == reference['n_neighbours']).all()
        
    def test_subset_lipids(self, universe):
        
        neighbours = Neighbours(universe, lipid_sel="name C", cutoff=10)
        neighbours.run()
    
        # it's a hexagonal lattice, but each hexagon is irregular (two sides longer than the other 4)
        # with cutoff=10, every lipid should have 2 neighbours
        reference = {
            'n_frames': 1,
            'n_residues': 50,
            'n_neighbours': np.array(
                [
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 0, 1, 0, 1
                ]
            )
        }
        
        assert neighbours.neighbours.shape[0] == (reference['n_frames'])
        assert neighbours.neighbours[0].shape == (reference['n_residues'], reference['n_residues'])
        assert_array_equal(np.sum(neighbours.neighbours[0].toarray(), axis=0), reference['n_neighbours'])


class TestNeighboursExceptions:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    def test_Exceptions(self, universe):
            
        match = "'cutoff' must be greater than 0"
        with pytest.raises(ValueError, match=match):
            Neighbours(
                universe=universe,
                lipid_sel="name L C",
                cutoff=-1
            )


class TestNeighboursCount:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C',
        'cutoff': 12.0
    }
    
    @pytest.fixture(scope='class')
    def neighbours(self, universe):
        neighbours = Neighbours(universe, **self.kwargs)
        neighbours.run()
        return neighbours
    
    # the sequence of n_CHOL_neighbours and n_LIPID_neighbours
    # arise because - even though the atoms are arranged on a
    # hexagonal lattice - each residue has two atoms
    @pytest.fixture(scope='class')
    def reference(self):
        
        reference = {
            'n_residues': 100,
            'n_frames': 1,
            'n_columns': 6,
            'n_neighbours': 6,
            'n_CHOL_neighbours': np.array(
                [
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3,
                    3, 2, 4, 2, 5, 1, 4, 2, 4, 3
                ]
            )
        }
        reference['n_LIPID_neighbours'] = 6 - reference['n_CHOL_neighbours']
        
        return reference
    
    def test_count_neighbours(self, neighbours, reference):
        
        counts = neighbours.count_neighbours()
        
        assert_array_equal(counts.shape, (reference['n_residues'] * reference['n_frames'], reference['n_columns']))
        assert_array_equal(counts.Frame, np.zeros(reference['n_residues']))
        assert (counts.Total == reference['n_neighbours']).all()
        assert_array_equal(counts.nCHOL, reference['n_CHOL_neighbours'])
        assert_array_equal(counts.nLIPI, reference['n_LIPID_neighbours'])
    
    def test_count_neighbours_count_by(self, neighbours, reference):
        
        # make every LIPI take the value 0
        # and every CHOL take the value 1
        count_by = np.zeros((reference['n_residues'], reference['n_frames']), dtype=np.int8)
        count_by[1::2] = 1
        
        counts = neighbours.count_neighbours(count_by=count_by)
        
        assert_array_equal(counts.shape, (reference['n_residues'] * reference['n_frames'], reference['n_columns']))
        assert_array_equal(counts.Frame, np.zeros(reference['n_residues']))
        assert (counts.Total == reference['n_neighbours']).all()
        assert_array_equal(counts.n0, reference['n_LIPID_neighbours'])
        assert_array_equal(counts.n1, reference['n_CHOL_neighbours'])
        
    def test_count_neighbours_count_by_offset(self, neighbours, reference):
        
        # make every LIPI take the value 0
        # and every CHOL take the value 1
        count_by = np.zeros((reference['n_residues'], reference['n_frames']), dtype=np.int8)
        count_by[1::2] = 1
        
        # check it doesn't matter what numbers are in 'count_by'
        # make LIPI = 100, CHOL = 110
        count_by = (count_by + 10) * 10
        
        counts = neighbours.count_neighbours(count_by=count_by)
        
        assert_array_equal(counts.shape, (reference['n_residues'] * reference['n_frames'], reference['n_columns']))
        assert_array_equal(counts.Frame, np.zeros(reference['n_residues']))
        assert (counts.Total == reference['n_neighbours']).all()
        assert_array_equal(counts.n100, reference['n_LIPID_neighbours'])
        assert_array_equal(counts.n110, reference['n_CHOL_neighbours'])

    def test_count_neighbours_count_by_labels(self, neighbours, reference):
        
        # make every LIPI take the value 0
        # and every CHOL take the value 1
        count_by = np.zeros((reference['n_residues'], reference['n_frames']), dtype=np.int8)
        count_by[1::2] = 1
        
        # liquid-disordered  = 0
        # liquid-ordered = 1
        count_by_labels = {"Ld": 0, "Lo": 1}
        
        counts = neighbours.count_neighbours(count_by=count_by, count_by_labels=count_by_labels)
        
        assert_array_equal(counts.shape, (reference['n_residues'] * reference['n_frames'], reference['n_columns']))
        assert_array_equal(counts.Frame, np.zeros(reference['n_residues']))
        assert (counts.Total == reference['n_neighbours']).all()
        assert_array_equal(counts.nLd, reference['n_LIPID_neighbours'])
        assert_array_equal(counts.nLo, reference['n_CHOL_neighbours'])

    def test_run_method_not_called(self, universe):
        
        neighbours = Neighbours(
            universe=universe,
            lipid_sel='name L C',
            cutoff=12.0
        )
        
        match = ".neighbours attribute is None: use .run\(\) before calling .count_neighbours\(\)"  # noqa:W605
        with pytest.raises(NoDataError, match=match):
            neighbours.count_neighbours()


class TestNeighboursEnrichment:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C',
        'cutoff': 10.0
    }
    
    @pytest.fixture(scope='class')
    def neighbours(self, universe):
        neighbours = Neighbours(universe, **self.kwargs)
        neighbours.run()
        return neighbours
        
    reference = {
        'n_species': 2,
        'columns': np.array(["Label", "Frame", "feCHOL", "feLIPI"]),
        'enrichment': np.array(
            [
                [0.4, 1.6],
                [1.6, 0.4]
            ])
    }
    
    def test_enrichment_cutoff(self, neighbours):
        
        counts, enrichment = neighbours.count_neighbours(return_enrichment=True)
        
        assert_array_equal((self.reference['n_species'], self.reference['columns'].size), enrichment.shape)
        assert_array_equal(self.reference['columns'], enrichment.columns)
        assert_array_equal(self.reference['enrichment'], enrichment[["feCHOL", "feLIPI"]].values)
    

class TestNeighboursClusters:
    
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(HEX_LAT)

    kwargs = {
        'lipid_sel': 'name L C',
        'cutoff': 12.0
    }
    
    @pytest.fixture(scope='class')
    def neighbours(self, universe):
        neighbours = Neighbours(universe, **self.kwargs)
        neighbours.run()
        return neighbours
    
    @pytest.fixture(scope='class')
    def reference(self):
        
        reference = {
            'n_frames': 1,
            'largest_leaflet': 50,
            'largest_lipid': 15,
            'largest_chol': 15,
            'all_upper_indices': np.arange(50),
            'all_lower_indices': np.arange(50, 100),
            'lipid_upper_indices': np.array(
                [
                    0, 4, 6, 10, 14, 16, 20, 24, 26, 30,
                    34, 36, 40, 44, 46
                ]
            ),
            'chol_lower_indices': np.array(
                [
                    53, 55, 59, 63, 65, 69, 73, 75, 79, 83,
                    85, 89, 93, 95, 99
                ]
            ),
            'lipid_sel': "name L",
            'chol_sel': 'name C',
            'all_sel': 'name L C',
            'leaflets': np.ones(100, dtype=np.int8)
        }
        
        # first 50 residues in the upper leaflet
        # next 50 in the lower leaflet
        reference['leaflets'][50:] = -1
        
        return reference
    
    def test_clusters(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel=reference['all_sel'],
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_leaflet'])
        assert_array_equal(largest_cluster_indices[0], reference['all_upper_indices'])
        
    def tests_clusters_lower_leaflet(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel=reference['all_sel'],
            filter_by=reference['leaflets'] == -1,
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_leaflet'])
        assert_array_equal(largest_cluster_indices[0], reference['all_lower_indices'])
        
    def test_clusters_upper_leaflet_lipid(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel=reference['lipid_sel'],
            filter_by=reference['leaflets'] == 1,
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_lipid'])
        assert_array_equal(largest_cluster_indices[0], reference['lipid_upper_indices'])
        
    def test_clusters_lower_leaflet_chol(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel=reference['chol_sel'],
            filter_by=reference['leaflets'] == -1,
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_chol'])
        assert_array_equal(largest_cluster_indices[0], reference['chol_lower_indices'])
        
    def test_clusters_dont_return_resindices(self, neighbours, reference):
        
        largest_cluster = neighbours.largest_cluster(
            cluster_sel=reference['all_sel'],
            filter_by=reference['leaflets'] == 1,
            return_indices=False
        )
        
        assert isinstance(largest_cluster, np.ndarray)
        assert largest_cluster.size == reference['n_frames']
        
    def test_no_cluster_sel(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            filter_by=reference['leaflets'] == 1,
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_leaflet'])
        assert_array_equal(largest_cluster_indices[0], reference['all_upper_indices'])
        
    def test_no_filter_by(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel="name L C",
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_leaflet'])
        assert_array_equal(largest_cluster_indices[0], reference['all_upper_indices'])
        
    def test_filter_by_2D(self, neighbours, reference):
        
        largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
            cluster_sel="name L C",
            filter_by=reference['leaflets'][:, np.newaxis] == 1,
            return_indices=True
        )
        
        assert largest_cluster.size == reference['n_frames']
        assert_array_equal(largest_cluster, reference['largest_leaflet'])
        assert_array_equal(largest_cluster_indices[0], reference['all_upper_indices'])
        
    def test_run_method_not_called(self, universe):
        
        neighbours = Neighbours(
            universe=universe,
            lipid_sel='name L C',
            cutoff=12.0
        )
        
        match = ".neighbours attribute is None: use .run\(\) before calling .largest_cluster\(\)"  # noqa:W605
        with pytest.raises(NoDataError, match=match):
            neighbours.largest_cluster()

    def test_bad_cluster_sel(self, neighbours, reference):
          
        match = "'cluster_sel' produces atom empty AtomGroup. Please check the selection string."
        with pytest.raises(ValueError, match=match):
            largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
                cluster_sel="",
                filter_by=reference['leaflets'] == -1,
                return_indices=True
            )
            
    def test_bad_filter_by_dimensions(self, neighbours, reference):
          
        match = "'filter_by' must either be a 1D array containing non-changing boolean"
        with pytest.raises(ValueError, match=match):
            largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
                cluster_sel=reference['all_sel'],
                filter_by=np.array(None),
                return_indices=True
            )
            
    def test_bad_filter_num_lipids(self, neighbours, reference):
          
        match = "The shape of 'filter_by' must be \(n_residues,\)"  # noqa:W605
        with pytest.raises(ValueError, match=match):
            largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
                cluster_sel=reference['all_sel'],
                filter_by=reference['leaflets'][:99] == 1,
                return_indices=True
            )
