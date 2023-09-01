from ._lipyferric import __version__

__version__ = version = __version__.replace("-", ".")  # Rust uses X.Y-dev0 but Python X.Y.dev0
__version_tuple__ = version_tuple = tuple(version.split("."))

from lipyphilic import (
    analysis,
    plotting,
    transformations,
)
from lipyphilic.leaflets.assign_leaflets import (
    AssignCurvedLeaflets,
    AssignLeaflets,
)

__all__ = [
    "AssignLeaflets",
    "AssignCurvedLeaflets",
    "analysis",
    "plotting",
    "transformations",
]
