import warnings

import lipyphilic as lpp

__all__ = []

_msg = (
    "The module `lipyphilic.lib.flip_flop` is deprecated and will be removed "
    "in a later version. Please use `lipyphilic.analysis.flip_flop` instead."
)
warnings.warn(
    _msg,
    DeprecationWarning,
    stacklevel=2,
)


class FlipFlop:
    """The class is deprecated. Please use `lipyphilic.analysis.FlipFlop instead`."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.lib.flip_flop.FlipFlop` is deprecated and will be removed "
            "in a later version. Please use `lipyphilic.analysis.FlipFlop` instead."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        return lpp.analysis.FlipFlop(*args, **kwargs)  # noqa:  PLE0101
