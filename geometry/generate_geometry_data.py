"""The entry point the README documents: the generator is build_geometry.py.

Kept as a name because scripts, macros and notes call it; it forwards to the
detector build, which defaults to the CMS sample.
"""

from build_geometry import main

if __name__ == "__main__":
    raise SystemExit(main())
