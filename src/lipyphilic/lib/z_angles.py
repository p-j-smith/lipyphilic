import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.z_angles` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.z_angles` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class ZAngles:
    """The class is deprecated. Please use `lipyphilic.analysis.ZAngles instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.z_angles.ZAngles` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.ZAngles` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.ZAngles(*args, **kwargs)  # noqa:  PLE0101
