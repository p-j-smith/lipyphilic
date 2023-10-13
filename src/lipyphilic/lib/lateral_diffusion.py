import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.lateral_diffusion` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.lateral_diffusion` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class MSD:
    """The class is deprecated. Please use `lipyphilic.analysis.MSD instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.lateral_diffusion.MSD` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.MSD` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.MSD(*args, **kwargs)  # noqa:  PLE0101
