import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.order_parameter` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.order_parameter` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class SCC:
    """The class is deprecated. Please use `lipyphilic.analysis.SCC instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.order_parameter.SCC` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.SCC` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.SCC(*args, **kwargs)  # noqa:  PLE0101
