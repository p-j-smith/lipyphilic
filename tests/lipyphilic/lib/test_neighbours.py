
import pytest
import numpy as np
import MDAnalysis

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
