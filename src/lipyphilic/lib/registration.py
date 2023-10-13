import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.registration` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.registration` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class Registration:
    """The class is deprecated. Please use `lipyphilic.analysis.Registration instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.registration.Registration` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.Registration` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.Registration(*args, **kwargs)  # noqa:  PLE0101
