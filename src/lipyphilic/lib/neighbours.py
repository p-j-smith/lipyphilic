import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.neighbours` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.neighbours` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class Neighbours:
    """The class is deprecated. Please use `lipyphilic.analysis.Neighbours instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.neighbours.Neighbours` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.Neighbours` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.Neighbours(*args, **kwargs)  # noqa:  PLE0101
