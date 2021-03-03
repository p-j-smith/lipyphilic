"""
Simple systems created for testing lipyphilic
=============================================

"""
__all__ = [
    "HEX_LAT",
    "HEX_LAT_BUMP",
    "HEX_LAT_BUMP_MID_MOL",
    "HEX_LAT_BUMP_MID_ATOM",
    "HEX_LAT_OVERLAP",
    "ONE_CHOL",
    "ONE_CHOL_TRAJ",
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

# Same atoms as HEX_LAT but now each molecule has four atoms - two in the upper leaflet
# (z = 60) and two in the lower leaflet (z = 40).
HEX_LAT_MONO = resource_filename(__name__,
                                 "pdbs/HexGrid-4AtomsPerLipid-Monolayer.pdb")

# A single coarse-grained cholesterol molecule
ONE_CHOL = resource_filename(__name__,
                             "pdbs/Chol-Flip-Flop.pdb")
# A trajectory of 25 frames of the above cholesterol molecule
ONE_CHOL_TRAJ = resource_filename(__name__,
                                  "xtcs/Chol-Flip-Flop.xtc")
