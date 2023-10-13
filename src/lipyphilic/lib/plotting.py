import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.plotting` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.plotting` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class ProjectionPlot:
    """The class is deprecated. Please use `lipyphilic.plotting.ProjectionPlot instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.plotting.ProjectionPlot` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.plotting.ProjectionPlot` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.plotting.ProjectionPlot(*args, **kwargs)  # noqa:  PLE0101


class JointDensity:
    """This class is deprecated. Please use `lipyphilic.plotting.JointDensity` instead.`"""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.plotting.JointDensity` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.plotting.JointDensity` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.plotting.JointDensity(*args, **kwargs)  # noqa:  PLE0101
