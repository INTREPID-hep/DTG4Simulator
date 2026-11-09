"""The detectors this build can write: one module each (see build_geometry)."""

from detectors import cms, prototype

__all__ = ["DETECTORS", "cms", "prototype"]

#: name -> module; add a file here and the CLI can build it
DETECTORS = {"cms": cms, "mutel": prototype}
