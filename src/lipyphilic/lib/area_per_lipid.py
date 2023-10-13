import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.area_per_lipid` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.area_per_lipid` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class AreaPerLipid:
    """This class is deprecated. Please use `lipyphilic.analysis.AreaPerLipid instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.area_per_lipid.AreaPerLipid` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.AreaPerLipid` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.AreaPerLipid(*args, **kwargs)  # noqa:  PLE0101
