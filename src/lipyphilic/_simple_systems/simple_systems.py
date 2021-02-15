"""
Simple systems created for testing lipyphilic
=============================================

"""
__all__ = [
    "HEX_LAT",
    "HEX_LAT_BUMP",
    "HEX_LAT_BUMP_MID_MOL",
    "HEX_LAT_BUMP_MID_ATOM",
    "HEX_LAT_OVERLAP"
]

from pkg_resources import resource_filename

# Two hexagonal lattices to represent upper and lower leaflets, separated in
# z by 20 Angstrom
# Alternating residues of LIPID and CHOL
# Two atoms per molecule
# Atoms are equally spaced
HEX_LAT = resource_filename(__name__,
                            "pdbs/HexGrid-2AtomsPerLipid.pdb")
# As immediately above but the central part of the lattices translated in
# z by +5 Angstrom
HEX_LAT_BUMP = resource_filename(__name__,
                                 "pdbs/HexGrid-2AtomsPerLipid-Undulating.pdb")
# As immediately above but one molecule (CHOL with resid 78) has been
# moved to the midplane
HEX_LAT_BUMP_MID_MOL = resource_filename(__name__,
                                         "pdbs/HexGrid-2AtomsPerLipid-Undulating-Midplane-Molecule.pdb")
# Same as HEX_LAT_BUMP but one atom of one molecule (CHOL with resid 78) has been
# moved to the midplane
HEX_LAT_BUMP_MID_ATOM = resource_filename(__name__,
                                          "pdbs/HexGrid-2AtomsPerLipid-Undulating-Midplane-Atom.pdb")
# Same as HEX_LAT but the first two atoms are overalpping (have to exact same position)
HEX_LAT_OVERLAP = resource_filename(__name__,
                                    "pdbs/HexGrid-2AtomsPerLipid-OverlappingAtoms.pdb")
