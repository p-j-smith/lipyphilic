
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
            'n_residues': 100,  # there are 200 atoms but 100 lipids in total
            'n_neighbours': 6,
        }
    
        assert neighbours.neighbours.shape == (reference['n_residues'], reference['n_residues'])
        assert (np.sum(neighbours.neighbours.toarray(), axis=0) == 6).all()
        
    def test_neighbours_cutoff10(self, universe):
        
        neighbours = Neighbours(universe, **self.kwargs, cutoff=10)
        neighbours.run()
    
        # it's a hexagonal lattice, but each hexagon is irregular (two sides longer than the other 4)
        # with cutoff=10, every lipid should have 2 neighbours
        reference = {
            'n_residues': 100,  # there are 200 atoms but 100 lipids in total
            'n_neighbours': 6,
        }
    
        assert neighbours.neighbours.shape == (reference['n_residues'], reference['n_residues'])
        assert (np.sum(neighbours.neighbours.toarray(), axis=0)).all()


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
            lipid_sel='name LC',
            cutoff=12.0
        )
        
        match = ".neighbours attribute is None: use .run\(\) before calling .count_neighbours\(\)"  # noqa:W605
        with pytest.raises(NoDataError, match=match):
            neighbours.count_neighbours()
