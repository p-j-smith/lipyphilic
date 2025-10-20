"""
Simple systems created for testing lipyphilic
=============================================

"""
from importlib.resources import as_file, files

__all__ = [
    "HEX_LAT",
    "HEX_LAT_BUMP",
    "HEX_LAT_BUMP_MID_ATOM",
    "HEX_LAT_BUMP_MID_MOL",
    "HEX_LAT_OVERLAP",
    "ONE_CHOL",
    "ONE_CHOL_TRAJ",
    "TRICLINIC",
]

pdb_dir = files("__name__").joinpath("pdbs")
xtc_dir = files("__name__").joinpath("xtcs")


# Two hexagonal lattices to represent upper and lower leaflets, separated in
# z by 20 Angstrom
# Alternating residues of LIPID and CHOL
# Two atoms per molecule
# Atoms are equally spaced
HEX_LAT = as_file(pdb_dir.joinpath("HexGrid-2AtomsPerLipid.pdb"))
# As immediately above but the central part of the lattices translated in
# z by +5 Angstrom
HEX_LAT_BUMP = as_file(pdb_dir.joinpath("HexGrid-2AtomsPerLipid-Undulating.pdb"))
# As immediately above but one molecule (CHOL with resid 78) has been
# moved to the midplane
HEX_LAT_BUMP_MID_MOL = as_file(pdb_dir.joinpath("HexGrid-2AtomsPerLipid-Undulating-Midplane-Molecule.pdb"))
# Same as HEX_LAT_BUMP but one atom of one molecule (CHOL with resid 78) has been
# moved to the midplane
HEX_LAT_BUMP_MID_ATOM = as_file(pdb_dir.joinpath("pdbs/HexGrid-2AtomsPerLipid-Undulating-Midplane-Atom.pdb"))
# Same as HEX_LAT but the first two atoms are overalpping (have to exact same position)
HEX_LAT_OVERLAP = as_file(pdb_dir.joinpath("HexGrid-2AtomsPerLipid-OverlappingAtoms.pdb"))
# Same as HEX_LAT but the bilayer is centered at z=0 and the atoms wrapped into the box,
# plus the leaflets are separated by 8 Angstrom in z rather than 20
HEX_LAT_SPLIT_Z = as_file(pdb_dir.joinpath("HexGrid-2AtomsPerLipid-SplitAcrossZ.pdb"))

# Same atoms as HEX_LAT but now each molecule has four atoms - two in the upper leaflet
# (z = 60) and two in the lower leaflet (z = 40).
HEX_LAT_MONO = as_file(pdb_dir.joinpath("HexGrid-4AtomsPerLipid-Monolayer.pdb"))

HEX_LAT_TRANS = as_file(pdb_dir.joinpath("pdbs/HexGrid-2AtomsPerLipid-TranslatedIn_y.pdb"))
HEX_LAT_TRANS_TRAJ = as_file(xtc_dir.joinpath("xtcs/HexGrid-2AtomsPerLipid-TranslatedIn_y.xtc"))

# A single coarse-grained cholesterol molecule
ONE_CHOL = as_file(pdb_dir.joinpath("Chol-Flip-Flop.pdb"))
# A trajectory of 25 frames of the above cholesterol molecule
ONE_CHOL_TRAJ = as_file(xtc_dir.joinpath("Chol-Flip-Flop.xtc"))

# Triclinic system with three atoms
# One atom is in the center of the box, the other two outside the primary unit cell
TRICLINIC = as_file(pdb_dir.joinpath("Triclinic-3Atoms.pdb"))
