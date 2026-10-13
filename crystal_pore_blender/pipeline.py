"""CIF ingestion, periodic geometry, and display construction.

All Cartesian coordinates use fractional @ a row-stored lattice in angstroms.
Display scaling never changes the numerical pore model.
"""

from pathlib import Path
import itertools
import json
import shutil
import subprocess
import warnings

import numpy as np
from pymatgen.core import Element, Lattice, Structure
from pymatgen.io.cif import CifParser
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

from .connectivity import periodic_components
from .constants import LEGACY_VDW


def dump(path, data):
    Path(path).write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def load_job(path):
    path = Path(path).resolve()
    job = json.loads(path.read_text(encoding="utf-8"))
    if "cif" not in job:
        raise ValueError("job.json needs a 'cif' path.")
    job["cif"] = str((path.parent / job["cif"]).resolve())
    if job.get("zeopp_executable") and not shutil.which(job["zeopp_executable"]):
        job["zeopp_executable"] = str((path.parent / job["zeopp_executable"]).resolve())
    return job


def intake(cif, A=None):
    # Disable fractional rounding: preserve the coordinates written by the author.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        parser = CifParser(cif, occupancy_tolerance=1.0, frac_tolerance=0)
        structures = parser.parse_structures(primitive=False)
    if len(structures) != 1:
        raise ValueError(
            "Supply a CIF containing exactly one structural model/data block."
        )
    source = structures[0]
    if any(
        not site.is_ordered or abs(sum(site.species.values()) - 1) > 1e-6
        for site in source
    ):
        raise ValueError(
            "Disorder/partial occupancy: supply an explicitly ordered model."
        )
    lattice = np.asarray(source.lattice.matrix if A is None else A, dtype=float)
    if (
        lattice.shape != (3, 3)
        or not np.isfinite(lattice).all()
        or np.linalg.det(lattice) <= 1e-8
    ):
        raise ValueError("Use a finite, right-handed 3 x 3 lattice.")
    if not np.allclose(
        lattice @ lattice.T,
        source.lattice.matrix @ source.lattice.matrix.T,
        atol=1e-4,
        rtol=1e-6,
    ):
        raise ValueError(
            "Display lattice changes cell metrics; only Cartesian reorientation is allowed."
        )
    atoms = [
        {
            "source_index": i,
            "element": site.specie.symbol,
            "label": site.label,
            "fractional": (site.frac_coords % 1).tolist(),
        }
        for i, site in enumerate(source)
    ]
    return source, {
        "source_cif": str(Path(cif).resolve()),
        "lattice_A": lattice.tolist(),
        "atoms": atoms,
        "ordered": True,
        "parser_warnings": [str(w.message) for w in caught],
    }


def radii(job, elements):
    model = job.get("radii_model", "pymatgen")
    if model not in ("pymatgen", "legacy", "custom"):
        raise ValueError("radii_model must be pymatgen, legacy, or custom.")
    values = {}
    for symbol in set(elements):
        value = (
            Element(symbol).van_der_waals_radius
            if model == "pymatgen"
            else LEGACY_VDW.get(symbol)
            if model == "legacy"
            else None
        )
        if value is not None:
            values[symbol] = float(value)
    values.update(job.get("vdw_radii_A", {}))
    missing = set(elements) - values.keys()
    if missing:
        raise ValueError(
            f"Set vdw_radii_A for elements without tabulated radii: {sorted(missing)}"
        )
    values = {e: float(values[e]) for e in set(elements)}
    if any(not np.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError("vdW radii must be finite and positive.")
    return values


def validate_job(job):
    source, data = intake(job["cif"], job.get("lattice_A"))
    groups = job.get("guest_groups", [])
    if any(not isinstance(g, list) or not g for g in groups):
        raise ValueError("Each guest group must be a nonempty list of source indices.")
    removed = [i for g in groups for i in g]
    if any(type(i) is not int or not 0 <= i < len(source) for i in removed) or len(
        set(removed)
    ) != len(removed):
        raise ValueError(
            "Invalid or overlapping guest indices; inspect this CIF first."
        )
    if len(removed) == len(source):
        raise ValueError("No framework remains.")
    r = radii(job, [a["element"] for a in data["atoms"]])
    for key in ("probe_radius_A", "grid_step_A", "convergence_step_A"):
        if key in job and (not np.isfinite(float(job[key])) or float(job[key]) <= 0):
            raise ValueError(f"{key} must be finite and positive.")
    if job.get("surface", "contact") not in ("contact", "probe-center"):
        raise ValueError("surface must be contact or probe-center.")
    if job.get("accessibility", "percolating") not in ("percolating", "all"):
        raise ValueError("accessibility must be percolating or all.")
    if job.get("display_mode", "cell") not in ("cell", "single"):
        raise ValueError("display_mode must be cell or single.")
    if job.get("framework_mode", "network") not in ("network", "molecules"):
        raise ValueError("framework_mode must be network or molecules.")
    if "crop_fractional" in job:
        lo, hi = crop_bounds(job)
    if "pore_axis_uvw" in job:
        axis = np.asarray(job["pore_axis_uvw"], float)
        if (
            axis.shape != (3,)
            or not np.isfinite(axis).all()
            or np.linalg.norm(axis) == 0
        ):
            raise ValueError("pore_axis_uvw needs a nonzero finite three-vector.")
    for record in job.get("guest_images", []):
        group = record.get("group")
        shift = np.asarray(record.get("shift", []))
        if (
            type(group) is not int
            or not 0 <= group < len(groups)
            or shift.shape != (3,)
            or not np.isfinite(shift.astype(float)).all()
            or not np.equal(shift, np.rint(shift)).all()
        ):
            raise ValueError(
                "guest_images needs an existing group and three integer translations."
            )
    return source, data, removed, r


def crop_bounds(job):
    crop = job.get("crop_fractional", {"min": [0, 0, 0], "max": [1, 1, 1]})
    lo, hi = np.asarray(crop["min"], float), np.asarray(crop["max"], float)
    if (
        lo.shape != (3,)
        or hi.shape != (3,)
        or not np.isfinite([lo, hi]).all()
        or not np.all(hi > lo)
    ):
        raise ValueError(
            "crop_fractional needs finite min/max three-vectors with max > min."
        )
    return lo, hi


def replicas(A, reach):
    heights = 1 / np.linalg.norm(np.linalg.inv(A), axis=0)
    pads = np.maximum(1, np.ceil(reach / heights).astype(int))
    return np.asarray(list(itertools.product(*[range(-n, n + 1) for n in pads])))


def query_chunks(tree, X, cutoff, chunk=100000):
    result = np.empty(len(X))
    for start in range(0, len(X), chunk):
        result[start : start + chunk] = tree.query(
            X[start : start + chunk], distance_upper_bound=cutoff
        )[0]
    return result


def fields(F, elements, A, r, probe, step, max_voxels, accessibility="percolating"):
    N = np.maximum(2, np.ceil(np.linalg.norm(A, axis=1) / step).astype(int))
    count = int(np.prod(N, dtype=np.int64))
    if count > max_voxels:
        raise ValueError(
            f"Grid {N.tolist()} ({count:,} voxels) exceeds max_voxels={max_voxels:,}; increase grid_step_A."
        )
    frac = np.indices(N).reshape(3, -1).T / N
    X = frac @ A
    # Bound the longest grid-cell diagonal, including tilted cells.
    margin = sum(np.linalg.norm(A, axis=1) / N)
    clearance = np.full(len(X), probe + margin)
    elements = np.asarray(elements)
    for e in sorted(set(elements)):
        reach = r[e] + probe + margin
        translations = replicas(A, reach)
        points = (F[elements == e][None, :, :] + translations[:, None, :]).reshape(
            -1, 3
        ) @ A
        clearance = np.minimum(
            clearance, query_chunks(cKDTree(points), X, reach) - r[e]
        )
    clearance = clearance.reshape(N).astype(np.float32)
    _, percolating, components = periodic_components(clearance >= probe)
    centers_mask = percolating if accessibility == "percolating" else clearance >= probe
    if not centers_mask.any():
        raise ValueError(
            "No selected probe-center space. Try a smaller probe, a prepared guest-free model, or accessibility='all' for isolated cavities."
        )
    # Keep one tree of sampled centers rather than 27+ copies of a large grid.
    center_tree = cKDTree(frac[centers_mask.ravel()] @ A)
    distance = np.full(len(X), probe + margin)
    for T in replicas(A, probe + margin):
        distance = np.minimum(
            distance, query_chunks(center_tree, X - T @ A, probe + margin)
        )
    contact = (probe - distance).reshape(N).astype(np.float32)
    # Mask out isolated centers when drawing a percolating probe-center surface.
    center_level = clearance - probe
    center_level[(center_level >= 0) & ~centers_mask] = -margin
    return clearance, centers_mask, contact, center_level, components


def compute(job, out):
    source, data, removed, all_radii = validate_job(job)
    A = np.asarray(data["lattice_A"])
    retained = [i for i in range(len(source)) if i not in set(removed)]
    F = np.asarray([data["atoms"][i]["fractional"] for i in retained])
    els = [data["atoms"][i]["element"] for i in retained]
    r = {e: all_radii[e] for e in set(els)}
    probe, step = (
        float(job.get("probe_radius_A", 1.2)),
        float(job.get("grid_step_A", 0.25)),
    )
    accessibility = job.get("accessibility", "percolating")
    _, mask, contact, center, components = fields(
        F, els, A, r, probe, step, int(job.get("max_voxels", 4000000)), accessibility
    )
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    report = {
        "passed": False,
        "probe_radius_A": probe,
        "grid_counts": list(map(int, contact.shape)),
        "max_grid_step_A": float(
            max(np.linalg.norm(A, axis=1) / np.asarray(contact.shape))
        ),
        "radii_model": job.get("radii_model", "pymatgen"),
        "vdw_radii_A": r,
        "source_atoms": len(source),
        "framework_atoms": len(retained),
        "removed_guest_indices": removed,
        "source_lattice_preserved": True,
        "accessibility": accessibility,
        "contact_fraction": float((contact > 0).mean()),
        "center_fraction": float(mask.mean()),
        "periodic_components": components,
        "clearance_note": "Capped clearance cannot measure maximum cavity diameter.",
        "zeopp": {"ran": False},
    }
    if "convergence_step_A" in job:
        _, coarse_mask, coarse, _, coarse_components = fields(
            F,
            els,
            A,
            r,
            probe,
            float(job["convergence_step_A"]),
            int(job.get("max_voxels", 4000000)),
            accessibility,
        )
        delta = abs(float((coarse > 0).mean()) - report["contact_fraction"])
        report["fraction_resolution_difference"] = delta
        report["coarse_center_fraction"] = float(coarse_mask.mean())
        report["coarse_periodic_components"] = coarse_components
        if delta > float(job.get("convergence_tolerance", 0.025)):
            dump(out / "calculation.json", report)
            raise ValueError("Contact volume has not converged; refine grid_step_A.")
    dump(out / "intake.json", data)
    Structure(
        source.lattice,
        [source[i].specie for i in retained],
        [source[i].frac_coords for i in retained],
    ).to(filename=str(out / "framework_for_pores.cif"))
    np.savez_compressed(
        out / "periodic_probe_field.npz",
        accessible_centers=mask,
        contact_level_A=contact,
        center_level_A=center,
        lattice_A=A,
        probe_radius_A=probe,
    )
    executable = job.get("zeopp_executable")
    if executable:
        # Always provide the exact numerical radii table to the independent tool.
        radii_path = (out / "zeopp_radii.rad").resolve()
        radii_path.write_text(
            "\n".join(f"{e} {v}" for e, v in sorted(r.items())) + "\n", encoding="utf-8"
        )
        executable = shutil.which(executable) or executable
        command = [
            executable,
            "-r",
            str(radii_path),
            "-ha",
            "-res",
            str((out / "zeopp.res").resolve()),
            "-chan",
            str(probe),
            str((out / "zeopp.chan").resolve()),
            str((out / "framework_for_pores.cif").resolve()),
        ]
        for stale in ("zeopp.res", "zeopp.chan"):
            (out / stale).unlink(missing_ok=True)
        result = subprocess.run(command, text=True, capture_output=True, timeout=300)
        (out / "zeopp.log").write_text(
            result.stdout + "\n" + result.stderr, encoding="utf-8"
        )
        good = result.returncode == 0 and all(
            (out / f).is_file() and (out / f).stat().st_size
            for f in ("zeopp.res", "zeopp.chan")
        )
        report["zeopp"] = {
            "ran": True,
            "passed": bool(good),
            "returncode": result.returncode,
            "command": command,
        }
        if not good:
            dump(out / "calculation.json", report)
            raise ValueError("Zeo++ failed; inspect zeopp.log.")
    report["passed"] = True
    dump(out / "calculation.json", report)
    return report


def gentle_smooth(V, F):
    edges = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    unique, count = np.unique(edges, axis=0, return_counts=True)
    boundary = np.unique(unique[count == 1])
    matrix = coo_matrix(
        (
            np.ones(len(unique) * 2),
            (np.r_[unique[:, 0], unique[:, 1]], np.r_[unique[:, 1], unique[:, 0]]),
        ),
        shape=(len(V), len(V)),
    ).tocsr()
    degree = np.maximum(1, np.asarray(matrix.sum(1)).ravel())
    X = V.copy()
    for _ in range(3):
        for strength in (0.25, -0.26):
            delta = matrix @ X / degree[:, None] - X
            delta[boundary] = 0
            X += strength * delta
    return X


def write_obj(path, V, F):
    with Path(path).open("w", encoding="utf-8") as stream:
        for p in V:
            stream.write("v " + " ".join(f"{v:.7f}" for v in p) + "\n")
        for triangle in F:
            stream.write("f " + " ".join(str(int(i) + 1) for i in triangle) + "\n")


def bond_limit(e1, e2, job):
    overrides = {
        tuple(sorted(k.split("-"))): float(v)
        for k, v in job.get("bond_cutoffs_A", {}).items()
    }
    pair = tuple(sorted((e1, e2)))
    if pair in overrides:
        return overrides[pair]
    nodes, donors = (
        set(job.get("node_elements", [])),
        set(job.get("donor_elements", ["N", "O"])),
    )
    if (e1 in nodes and e2 in donors) or (e2 in nodes and e1 in donors):
        return float(job.get("coordination_cutoff_A", 2.6))
    from pymatgen.analysis.molecule_structure_comparator import CovalentRadius

    a, b = CovalentRadius.radius.get(e1), CovalentRadius.radius.get(e2)
    return 1.15 * (a + b) if a is not None and b is not None else 0


def unwrap_group(ids, atoms, lattice, job):
    """Unwrap a finite molecule by its periodic bond graph, not nearest anchor images."""
    F = np.asarray([atoms[i]["fractional"] for i in ids])
    edges = {i: [] for i in range(len(ids))}
    for i, j in itertools.combinations(range(len(ids)), 2):
        distance, T = lattice.get_distance_and_image(F[i], F[j])
        if (
            0.35
            < distance
            <= bond_limit(atoms[ids[i]]["element"], atoms[ids[j]]["element"], job)
        ):
            edges[i].append((j, T))
            edges[j].append((i, -T))
    shifts, stack = {0: np.zeros(3, int)}, [0]
    while stack:
        i = stack.pop()
        for j, T in edges[i]:
            proposed = shifts[i] + T
            if j in shifts and not np.array_equal(shifts[j], proposed):
                raise ValueError(
                    "Guest group is a periodic network, not a finite molecule."
                )
            if j not in shifts:
                shifts[j] = proposed
                stack.append(j)
    if len(shifts) != len(ids):
        raise ValueError(
            "Guest group is disconnected with current bond cutoffs; adjust bond_cutoffs_A or split groups."
        )
    return F + np.asarray([shifts[i] for i in range(len(ids))])


def molecular_templates(atoms, lattice, job, excluded):
    ids = [
        i
        for i, a in enumerate(atoms)
        if i not in excluded and a["element"] not in job.get("node_elements", [])
    ]
    adjacency = {i: [] for i in ids}
    for i, j in itertools.combinations(ids, 2):
        d, T = lattice.get_distance_and_image(
            atoms[i]["fractional"], atoms[j]["fractional"]
        )
        if 0.35 < d <= bond_limit(atoms[i]["element"], atoms[j]["element"], job):
            adjacency[i].append(j)
            adjacency[j].append(i)
    seen, result = set(), []
    for start in ids:
        if start in seen:
            continue
        group, stack = [], [start]
        while stack:
            i = stack.pop()
            if i in seen:
                continue
            seen.add(i)
            group.append(i)
            stack.extend(adjacency[i])
        result.append((group, unwrap_group(group, atoms, lattice, job)))
    return result


def display(job, calculation, out):
    from skimage.measure import marching_cubes

    _, current, removed, r = validate_job(job)
    calculation, out = Path(calculation), Path(out)
    data = json.loads((calculation / "intake.json").read_text(encoding="utf-8"))
    report = json.loads((calculation / "calculation.json").read_text(encoding="utf-8"))
    if not report.get("passed"):
        raise ValueError("Calculation did not pass.")
    if current["atoms"] != data["atoms"] or not np.allclose(
        current["lattice_A"], data["lattice_A"], atol=1e-5
    ):
        raise ValueError("CIF/coordinates do not match cached field intake.")
    if (
        sorted(removed) != sorted(report["removed_guest_indices"])
        or any(abs(r[e] - v) > 1e-8 for e, v in report["vdw_radii_A"].items())
        or abs(float(job.get("probe_radius_A", 1.2)) - report["probe_radius_A"]) > 1e-8
        or job.get("accessibility", "percolating") != report["accessibility"]
    ):
        raise ValueError(
            "Probe/radii/guest selection/accessibility differs from cached calculation."
        )
    with np.load(calculation / "periodic_probe_field.npz") as z:
        A = z["lattice_A"]
        surface = job.get("surface", "contact")
        C = z["contact_level_A" if surface == "contact" else "center_level_A"]
    atoms, lattice = data["atoms"], Lattice(A)
    lo, hi = crop_bounds(job)
    N = np.asarray(C.shape)
    start, stop = np.floor(lo * N).astype(int), np.ceil(hi * N).astype(int)
    shape = stop - start + 1
    if int(np.prod(shape)) > int(job.get("max_display_voxels", 8000000)):
        raise ValueError(
            "Display crop exceeds max_display_voxels; reduce crop or resolution."
        )
    axes = [np.arange(a, b + 1) for a, b in zip(start, stop)]
    field = C[np.ix_(*[x % N[i] for i, x in enumerate(axes)])].copy()
    guests = []
    images = job.get(
        "guest_images",
        [
            {"group": g, "shift": [0, 0, 0]}
            for g in range(len(job.get("guest_groups", [])))
        ],
    )
    for record in images:
        ids = job["guest_groups"][record["group"]]
        base = unwrap_group(ids, atoms, lattice, job) + np.asarray(record["shift"])
        molecular = [
            {
                "source_index": i,
                "element": atoms[i]["element"],
                "fractional": f.tolist(),
                "position_A": (f @ A).tolist(),
            }
            for i, f in zip(ids, base)
        ]
        bonds = make_bonds(molecular, job)
        guests.append(
            {
                "source_guest": record["group"],
                "periodic_translation": record["shift"],
                "atoms": molecular,
                "bonds": bonds,
            }
        )
    labels, n = label(field > 0)
    mode = job.get("display_mode", "cell")
    anchors = np.asarray([g["atoms"][0]["fractional"] for g in guests])
    if mode == "single":
        if not len(anchors):
            if "pore_anchor_fractional" not in job:
                raise ValueError(
                    "Single-pore mode needs pore_anchor_fractional or guest_images."
                )
            anchors = np.asarray([job["pore_anchor_fractional"]], float)
        indices = np.rint(anchors * N - start).astype(int)
        if np.any(indices < 0) or np.any(indices >= shape):
            raise ValueError("Guest/pore anchor outside crop.")
        keep = [int(labels[tuple(i)]) for i in indices]
        if 0 in keep or len(set(keep)) != 1:
            raise ValueError(
                "Anchors are not inside ONE pore component; select joint periodic images/crop."
            )
        selected = labels == keep[0]
        for d in set(range(3)) - set(job.get("open_crop_axes", [0])):
            if any(np.take(selected, k, axis=d).any() for k in (0, -1)):
                raise ValueError("Crop truncates lateral pore wall; expand the bounds.")
        field[(field > 0) & ~selected] = -float(job.get("grid_step_A", 0.25))
    if not (field.min() < 0 < field.max()):
        raise ValueError(
            "No surface crosses this crop; change probe/crop/grid resolution."
        )
    V, F, _, _ = marching_cubes(field, level=0, allow_degenerate=False)
    V = ((V + start) / N) @ A
    # Marching cubes default normals point toward positive values, i.e. the pore
    # interior. Reverse for an outward-facing void boundary.
    F = F[:, ::-1]
    smooth = gentle_smooth(V, F) if job.get("smooth_display", True) else V.copy()
    motion = np.linalg.norm(smooth - V, axis=1)
    if motion.max() > float(job.get("max_display_smoothing_A", 0.15)):
        raise ValueError("Display smoothing exceeds allowed motion.")
    records, excluded = [], set(removed)
    cutoff = float(job.get("wall_extra_A", 0.45))
    tree = cKDTree(V)
    inv = np.linalg.inv(A)
    templates = (
        molecular_templates(atoms, lattice, job, excluded)
        if job.get("framework_mode", "network") == "molecules"
        else [
            ([i], np.asarray([a["fractional"]]))
            for i, a in enumerate(atoms)
            if i not in excluded
        ]
    )
    # In molecule mode, coordination nodes are separate periodic records.
    if job.get("framework_mode", "network") == "molecules":
        templates += [
            ([i], np.asarray([a["fractional"]]))
            for i, a in enumerate(atoms)
            if i not in excluded and a["element"] in job.get("node_elements", [])
        ]
    for ids, base in templates:
        padding = (max(r.values()) + cutoff) * np.linalg.norm(inv, axis=0)
        lower = np.floor(lo - base.max(0) - padding).astype(int)
        upper = np.ceil(hi - base.min(0) + padding).astype(int)
        if int(np.prod(upper - lower + 1)) > int(
            job.get("max_replicas_per_group", 10000)
        ):
            raise ValueError(
                "Too many periodic display replicas; reduce crop or use a better conditioned cell."
            )
        for T in itertools.product(*[range(a, b + 1) for a, b in zip(lower, upper)]):
            X = (base + T) @ A
            if not np.any(
                tree.query(X)[0]
                < np.asarray([r[atoms[i]["element"]] + cutoff for i in ids])
            ):
                continue
            for i, f, x in zip(ids, base + T, X):
                records.append(
                    {
                        "source_index": i,
                        "element": atoms[i]["element"],
                        "position_A": x.tolist(),
                        "fractional": f.tolist(),
                        "role": "framework",
                    }
                )
            if len(records) > int(job.get("max_display_atoms", 30000)):
                raise ValueError("Too many display atoms; reduce crop or wall_extra_A.")
    if records and cKDTree([a["position_A"] for a in records]).query_pairs(0.2):
        raise ValueError("Duplicate/overlapping displayed framework sites.")
    bonds = make_bonds(records, job)
    coordinate_error = 0.0
    for atom in records + [a for g in guests for a in g["atoms"]]:
        diff = (
            np.asarray(atom["fractional"]) - atoms[atom["source_index"]]["fractional"]
        )
        coordinate_error = max(
            coordinate_error, float(np.linalg.norm((diff - np.rint(diff)) @ A))
        )
    if coordinate_error > 1e-5:
        raise ValueError(
            "Display coordinates changed beyond integer lattice translations."
        )
    axis = job.get("pore_axis_uvw")
    percolating = [c for c in report["periodic_components"] if c["percolates"]]
    # An inferred axis is meaningful only for a single rank-one channel family.
    if axis is None and len(percolating) == 1 and percolating[0]["dimension"] == 1:
        axis = percolating[0]["winding_vectors"][0]
    if not len(anchors):
        anchors = np.asarray([(lo + hi) / 2])
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out / "pore_mesh.npz", raw_vertices_A=V, vertices_A=smooth, triangles=F
    )
    write_obj(out / "pore_raw.obj", V, F)
    write_obj(out / "pore_display.obj", smooth, F)
    dump(
        out / "display.json",
        {
            "framework_atoms": records,
            "framework_bonds": bonds,
            "guests": guests,
            "lattice_A": A.tolist(),
            "pore_axis_uvw": axis,
            "probe_radius_A": report["probe_radius_A"],
            "surface": surface,
            "vdw_radii_A": r,
            "guest_anchor_position_A": (anchors @ A).mean(0).tolist(),
            "pore_components": 1 if mode == "single" else n,
            "display_mode": mode,
        },
    )
    validation = {
        "passed": True,
        "framework_atoms": len(records),
        "guest_molecules": len(guests),
        "display_mode": mode,
        "coordinate_modulo_lattice_error_A": coordinate_error,
        "smoothing_max_motion_A": float(motion.max()),
        "smoothing_mean_motion_A": float(motion.mean()),
        "triangles": len(F),
        "surface": surface,
    }
    dump(out / "validation.json", validation)
    return validation


def make_bonds(records, job):
    if not records:
        return []
    elements = sorted(set(a["element"] for a in records))
    limits = [bond_limit(a, b, job) for a in elements for b in elements]
    X = np.asarray([a["position_A"] for a in records])
    result = []
    for i, j in sorted(cKDTree(X).query_pairs(max(limits, default=0))):
        distance = float(np.linalg.norm(X[i] - X[j]))
        if (
            0.35
            < distance
            <= bond_limit(records[i]["element"], records[j]["element"], job)
        ):
            result.append([i, j, distance])
    return result
