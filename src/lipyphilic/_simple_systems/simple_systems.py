"""
Simple systems created for testing lipyphilic
=============================================

"""
__all__ = [
    "HEX_LAT",
    "HEX_LAT_BUMP"
]

from pkg_resources import resource_filename

HEX_LAT = resource_filename(__name__,
                            "pdbs/HexGrid-2AtomsPerLipid.pdb")
HEX_LAT_BUMP = resource_filename(__name__,
                                 "pdbs/HexGrid-2AtomsPerLipid-Undulating.pdb")
