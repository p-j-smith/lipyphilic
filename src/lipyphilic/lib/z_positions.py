import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.z_positions` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.z_positions` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class ZPositions:
    """The class is deprecated. Please use `lipyphilic.analysis.ZPositions instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.z_positions.ZPositions` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.ZPositions` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.ZPositions(*args, **kwargs)  # noqa:  PLE0101
