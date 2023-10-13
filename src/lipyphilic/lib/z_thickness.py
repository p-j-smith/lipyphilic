import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.z_thickness` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.z_thickness` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class ZThickness:
    """The class is deprecated. Please use `lipyphilic.analysis.ZThickness instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.z_thickness.ZThickness` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.ZThickness` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.ZThickness(*args, **kwargs)  # noqa:  PLE0101
