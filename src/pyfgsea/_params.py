"""Helpers for config-driven parameter dictionaries."""

from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Any, Optional


def merge_params(
    params: Optional[Mapping[str, Any]] = None,
    *,
    allowed_keys: Optional[Collection[str]] = None,
    context: str = "parameters",
    **kwargs: Any,
) -> dict[str, Any]:
    """Merge a params mapping with explicit keyword arguments.

    Explicit keyword arguments intentionally override values from ``params``.
    Unknown keys are rejected early when ``allowed_keys`` is supplied.
    """

    if params is None:
        merged: dict[str, Any] = {}
    elif isinstance(params, Mapping):
        merged = dict(params)
    else:
        raise TypeError(f"{context} params must be a mapping or None.")

    merged.update(kwargs)

    if allowed_keys is not None:
        unknown = sorted(set(merged) - set(allowed_keys))
        if unknown:
            allowed = ", ".join(sorted(allowed_keys))
            unknown_display = ", ".join(unknown)
            raise ValueError(f"Unknown {context} parameter(s): {unknown_display}. Allowed keys: {allowed}.")

    return merged


def pop_required(merged: dict[str, Any], name: str, value: Any, context: str) -> Any:
    """Resolve a required argument from a positional value or params mapping."""

    from_params = merged.pop(name, None)
    resolved = value if value is not None else from_params
    if resolved is None:
        raise TypeError(f"{context} requires '{name}' as an argument or params['{name}'].")
    return resolved
