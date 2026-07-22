"""Validated structure preparation for bilayer stacking scans."""

from .bilayer import (
    MANIFEST_SCHEMA,
    BilayerError,
    build_bilayer,
    generate_scan,
    periodic_shifts,
    read_monolayer,
    write_bilayer,
)

__all__ = [
    "MANIFEST_SCHEMA",
    "BilayerError",
    "build_bilayer",
    "generate_scan",
    "periodic_shifts",
    "read_monolayer",
    "write_bilayer",
]

__version__ = "0.1.0a1"
