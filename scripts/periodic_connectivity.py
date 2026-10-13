"""Compatibility import; implementation lives in the installed package."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from crystal_pore_blender.connectivity import periodic_components

__all__ = ["periodic_components"]
