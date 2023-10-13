import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.assign_leaflets` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.leaflets.assign_leaflets` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class AssignLeaflets:
    """The class is deprecated. Please use `lipyphilic.AssignLeaflets instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.assign_leaflets.AssignLeaflets` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.AssignLeaflets` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.AssignLeaflets(*args, **kwargs)  # noqa:  PLE0101


class AssignCurvedLeaflets:
    """This class is deprecated. Please use `lipyphilic.AssignCurvedLeaflets instead.`"""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.assign_leaflets.AssignCurvedLeaflets` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.AssignCurvedLeaflets` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.AssignCurvedLeaflets(*args, **kwargs)  # noqa:  PLE0101
