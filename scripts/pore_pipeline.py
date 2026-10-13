"""Compatibility entry point; prefer the installed crystal-pore-blender command."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from crystal_pore_blender.cli import main

if __name__ == "__main__":
    sys.exit(main())
