try:
    from lipyphilic._version import version as __version__
except ImportError:
    __version__ = "not-installed"

from lipyphilic.leaflets.assign_leaflets import (
    AssignCurvedLeaflets,
    AssignLeaflets,
)
