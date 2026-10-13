"""Compatibility entry point for Blender's --python option."""

from pathlib import Path
import runpy

runpy.run_path(
    str(Path(__file__).resolve().parents[1] / "crystal_pore_blender" / "render.py"),
    run_name="__main__",
)
