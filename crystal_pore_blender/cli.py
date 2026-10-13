"""Installed command-line interface. Heavy geometry imports are lazy."""

from pathlib import Path
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys

from . import __version__


def find_blender(explicit=None):
    if explicit:
        result = shutil.which(explicit) or str(Path(explicit).expanduser())
        if not Path(result).is_file():
            raise ValueError(f"Blender executable not found: {explicit}")
        return result
    if shutil.which("blender"):
        return shutil.which("blender")
    candidates = [Path("/Applications/Blender.app/Contents/MacOS/Blender")]
    for root in [os.environ.get("ProgramFiles"), os.environ.get("LOCALAPPDATA")]:
        if root:
            candidates += sorted(
                Path(root).glob("Blender Foundation/Blender*/blender.exe"), reverse=True
            )
    return next((str(p) for p in candidates if p.is_file()), None)


def doctor(explicit=None):
    dependencies = {
        name: importlib.util.find_spec(name) is not None
        for name in ("numpy", "scipy", "pymatgen", "skimage", "PIL")
    }
    result = {
        "version": __version__,
        "python": sys.version.split()[0],
        "geometry_dependencies": dependencies,
        "blender_executable": find_blender(explicit),
        "bpy_available": importlib.util.find_spec("bpy") is not None,
        "zeopp_executable": shutil.which("network"),
    }
    result["geometry_ready"] = all(dependencies.values())
    result["render_ready"] = bool(
        result["blender_executable"] or result["bpy_available"]
    )
    print(json.dumps(result, indent=2))
    return result


def initialize(cif, target, force=False):
    from .pipeline import dump, intake

    _, data = intake(cif)
    target = Path(target).resolve()
    if target.exists() and not force:
        raise ValueError(f"{target} exists; use --force to replace it.")
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        cif_path = os.path.relpath(Path(cif).resolve(), target.parent)
    except ValueError:  # Windows drives cannot be relativized.
        cif_path = str(Path(cif).resolve())
    job = {
        "cif": cif_path,
        "guest_groups": [],
        "node_elements": [],
        "radii_model": "pymatgen",
        "probe_radius_A": 1.2,
        "grid_step_A": 0.25,
        "surface": "contact",
        "accessibility": "percolating",
        "display_mode": "cell",
        "framework_mode": "network",
        "crop_fractional": {"min": [0, 0, 0], "max": [1, 1, 1]},
        "smooth_display": True,
    }
    dump(target, job)
    dump(target.with_name(target.stem + ".intake.json"), data)
    print(f"Created {target} ({len(data['atoms'])} indexed atoms).")
    return job


def make_demo(out):
    """Create a synthetic channel and a reproducible job without fetching data."""
    import numpy as np
    from pymatgen.core import Lattice, Structure
    from .pipeline import dump

    inputs = out / "input"
    inputs.mkdir(parents=True, exist_ok=True)
    coordinates = [
        [x, 6 + 3.7 * np.cos(a), 6 + 3.7 * np.sin(a)]
        for x in np.arange(0, 12, 1.5)
        for a in np.arange(0, 2 * np.pi, 2 * np.pi / 16)
    ]
    cif = inputs / "demo-channel.cif"
    Structure(
        Lattice.cubic(12),
        ["C"] * len(coordinates),
        coordinates,
        coords_are_cartesian=True,
    ).to(filename=str(cif))
    job = {
        "cif": str(cif),
        "probe_radius_A": 1.2,
        "grid_step_A": 0.35,
        "pore_axis_uvw": [1, 0, 0],
        "display_mode": "single",
        "framework_mode": "network",
        "pore_anchor_fractional": [0.5, 0.5, 0.5],
        "open_crop_axes": [0],
        "crop_fractional": {"min": [0, 0.27, 0.27], "max": [1, 0.73, 0.73]},
    }
    dump(inputs / "job.json", job)
    return job


def render(display_dir, out, args):
    # Check availability before expensive geometry in run().
    executable = find_blender(args.blender)
    script = Path(__file__).with_name("render.py")
    command = (
        [
            executable,
            "--background",
            "--factory-startup",
            "--disable-autoexec",
            "--python-exit-code",
            "1",
            "--python",
            str(script),
            "--",
        ]
        if executable
        else [sys.executable, str(script), "--"]
    )
    command += [
        "--data",
        str(display_dir / "display.json"),
        "--mesh",
        str(display_dir / "pore_mesh.npz"),
        "--style",
        str(args.style or Path(__file__).parent / "styles" / "default.json"),
        "--out",
        str(out),
        "--mode",
        args.mode,
        "--view",
        args.view,
        "--section",
        args.section,
    ]
    subprocess.run(command, check=True)
    from .imaging import composite_white

    composite_white(out / "figure_transparent.png")


def run(job, out, args):
    from .pipeline import compute, display, dump

    if not args.no_render and not (
        find_blender(args.blender) or importlib.util.find_spec("bpy")
    ):
        raise ValueError(
            "Rendering needs Blender. Install Blender and add it to PATH, use --blender PATH, or install the optional [render] extra. Use --no-render to export geometry only."
        )
    if args.style and not Path(args.style).is_file():
        raise ValueError(f"Style file not found: {args.style}")
    if out.exists() and any(out.iterdir()) and not args.force:
        raise ValueError(
            f"Output directory {out} is not empty; choose another --out or use --force."
        )
    out.mkdir(parents=True, exist_ok=True)
    dump(out / "job_used.json", job)
    calculation = compute(job, out / "calculation")
    geometry = display(job, out / "calculation", out / "display")
    if not args.no_render:
        render(out / "display", out / "render", args)
    result = {
        "output": str(out),
        "calculation_passed": calculation["passed"],
        "display_passed": geometry["passed"],
        "rendered": not args.no_render,
    }
    dump(out / "run.json", result)
    print(json.dumps(result))
    return result


def parser():
    p = argparse.ArgumentParser(
        prog="crystal-pore-blender",
        description="Turn CIF crystal pores into editable Blender scenes.",
    )
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("doctor", help="Check Python dependencies and Blender discovery")
    d.add_argument("--blender")
    init = sub.add_parser("init", help="Create a generic job and indexed atom list")
    init.add_argument("--cif", required=True)
    init.add_argument("--out", type=Path, default=Path("job.json"))
    init.add_argument("--force", action="store_true")
    inspect = sub.add_parser("inspect", help="Expand symmetry and index atoms")
    inspect.add_argument("--cif", required=True)
    inspect.add_argument("--out", type=Path, required=True)
    for name in ("validate", "compute", "display", "run", "batch", "demo"):
        item = sub.add_parser(name)
        item.add_argument("--job", required=name not in ("run", "batch", "demo"))
        if name == "validate":
            continue
        item.add_argument("--out", type=Path, required=True)
        if name == "display":
            item.add_argument("--calculation", type=Path, required=True)
        if name in ("run", "batch", "demo"):
            item.add_argument(
                "--cif",
                required=name == "batch",
                nargs="+" if name == "batch" else None,
            )
            item.add_argument("--no-render", action="store_true")
            item.add_argument("--force", action="store_true")
            item.add_argument("--blender")
            item.add_argument("--style", type=Path)
            item.add_argument("--mode", choices=("full", "cutaway"), default="full")
            item.add_argument("--view", choices=("free", "side", "end"), default="free")
            item.add_argument(
                "--section",
                choices=("longitudinal", "transverse"),
                default="longitudinal",
            )
            item.add_argument("--probe", type=float)
            item.add_argument("--step", type=float)
            item.add_argument("--surface", choices=("contact", "probe-center"))
            item.add_argument("--accessibility", choices=("percolating", "all"))
            item.add_argument("--axis", type=float, nargs=3)
    item = sub.add_parser(
        "render", help="Render an existing display without recalculating"
    )
    item.add_argument("--display", type=Path, required=True)
    item.add_argument("--out", type=Path, required=True)
    item.add_argument("--blender")
    item.add_argument("--style", type=Path)
    item.add_argument("--mode", choices=("full", "cutaway"), default="full")
    item.add_argument("--view", choices=("free", "side", "end"), default="free")
    item.add_argument(
        "--section", choices=("longitudinal", "transverse"), default="longitudinal"
    )
    return p


def main():
    p = parser()
    args = p.parse_args()
    try:
        if args.command == "doctor":
            result = doctor(args.blender)
            return 0 if result["geometry_ready"] else 1
        from .pipeline import compute, display, dump, intake, load_job, validate_job

        if args.command == "init":
            initialize(args.cif, args.out, args.force)
        elif args.command == "inspect":
            _, data = intake(args.cif)
            args.out.mkdir(parents=True, exist_ok=True)
            dump(args.out / "intake.json", data)
            print(f"Indexed atoms: {len(data['atoms'])}")
        elif args.command == "render":
            if not (find_blender(args.blender) or importlib.util.find_spec("bpy")):
                raise ValueError(
                    "Blender is unavailable; run doctor or pass --blender PATH."
                )
            render(args.display.resolve(), args.out.resolve(), args)
        else:
            job = load_job(args.job) if args.job else {}
            if args.command in ("run", "batch", "demo"):
                if args.command == "run" and bool(args.cif) == bool(args.job):
                    raise ValueError("Use exactly one of --cif or --job.")
                for attr, key in (
                    ("probe", "probe_radius_A"),
                    ("step", "grid_step_A"),
                    ("surface", "surface"),
                    ("accessibility", "accessibility"),
                    ("axis", "pore_axis_uvw"),
                ):
                    if getattr(args, attr) is not None:
                        job[key] = getattr(args, attr)
                if args.command == "demo":
                    demo = make_demo(args.out.resolve())
                    demo.update(job)
                    run(demo, args.out.resolve() / "scene", args)
                elif args.command == "run":
                    if args.cif:
                        job["cif"] = str(Path(args.cif).resolve())
                    run(job, args.out.resolve(), args)
                else:
                    # Material-specific indices/crops cannot safely be reused across CIFs.
                    if job.get("guest_groups") or job.get("guest_images"):
                        raise ValueError(
                            "Batch templates cannot contain source-specific guest indices. Use individual jobs for guest selections."
                        )
                    records, seen, failed = [], set(), False
                    for cif in args.cif:
                        name = Path(cif).stem
                        if name in seen:
                            raise ValueError("Batch CIF basenames must be unique.")
                        seen.add(name)
                    args.out.mkdir(parents=True, exist_ok=True)
                    for cif in args.cif:
                        item = dict(job, cif=str(Path(cif).resolve()))
                        try:
                            records.append(
                                run(item, args.out.resolve() / Path(cif).stem, args)
                            )
                        except (
                            ValueError,
                            OSError,
                            subprocess.SubprocessError,
                        ) as error:
                            failed = True
                            records.append(
                                {"cif": cif, "passed": False, "error": str(error)}
                            )
                    dump(args.out / "batch.json", records)
                    if failed:
                        return 1
            elif args.command == "validate":
                _, data, removed, radii = validate_job(job)
                print(
                    json.dumps(
                        {
                            "valid": True,
                            "atoms": len(data["atoms"]),
                            "removed_atoms": len(removed),
                            "vdw_radii_A": radii,
                        }
                    )
                )
            elif args.command == "compute":
                print(json.dumps(compute(job, args.out.resolve())))
            else:
                print(
                    json.dumps(
                        display(job, args.calculation.resolve(), args.out.resolve())
                    )
                )
    except (
        ValueError,
        OSError,
        KeyError,
        ImportError,
        subprocess.SubprocessError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
