import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.memb_thickness` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.memb_thickness` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class MembThickness:
    """The class is deprecated. Please use `lipyphilic.analysis.MembThickness instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.memb_thickness.MembThickness` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.MembThickness` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.MembThickness(*args, **kwargs)  # noqa:  PLE0101
