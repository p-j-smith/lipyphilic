try:
    from lipyphilic._version import version as __version__
except ImportError:
    __version__ = "not-installed"

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
