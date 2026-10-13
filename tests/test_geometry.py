import json
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from crystal_pore_blender.connectivity import periodic_components
from crystal_pore_blender.pipeline import (
    compute,
    display,
    fields,
    intake,
    radii,
    unwrap_group,
)


def test_periodic_winding_distinguishes_cavity_and_channels():
    isolated = np.zeros((7, 7, 7), bool)
    isolated[2:5, 2:5, 2:5] = True
    _, mask, components = periodic_components(isolated)
    assert not mask.any()
    assert components[0]["dimension"] == 0
    channel = np.zeros_like(isolated)
    channel[:, 3, 3] = True
    _, mask, components = periodic_components(channel)
    assert mask.sum() == 7
    assert components[0]["dimension"] == 1
    assert components[0]["independent_pure_lattice_axes"] == [0]
    sheet = np.zeros_like(isolated)
    sheet[:, :, 3] = True
    assert periodic_components(sheet)[2][0]["dimension"] == 2
    assert periodic_components(np.ones_like(isolated))[2][0]["dimension"] == 3


def test_diagonal_winding_is_not_a_pure_axis():
    mask = np.zeros((7, 7, 3), bool)
    for i in range(7):
        mask[i, i, 1] = True
        mask[(i + 1) % 7, i, 1] = True
    component = periodic_components(mask)[2][0]
    assert component["dimension"] == 1
    assert component["independent_pure_lattice_axes"] == []
    assert abs(component["winding_vectors"][0][0]) == abs(
        component["winding_vectors"][0][1]
    )


def test_tilted_cell_clearance_matches_brute_periodic_distances():
    A = np.array([[6.0, 0, 0], [4.0, 5.0, 0], [1.0, 1.0, 6.0]])
    F = np.array([[0.2, 0.3, 0.4]])
    clearance, mask, contact, center, _ = fields(
        F, ["C"], A, {"C": 1.7}, 0.7, 1.0, 10000, "all"
    )
    N = np.array(clearance.shape)
    X = (np.indices(N).reshape(3, -1).T / N) @ A
    T = np.array(
        [(a, b, c) for a in range(-3, 4) for b in range(-3, 4) for c in range(-3, 4)]
    )
    points = (F + T) @ A
    exact = np.linalg.norm(X[:, None, :] - points[None, :, :], axis=2).min(1) - 1.7
    near = exact <= 0.7
    assert np.allclose(clearance.ravel()[near], exact[near], atol=1e-6)
    assert np.array_equal(mask.ravel(), exact >= 0.7)
    assert np.all(contact[mask] > 0)
    assert np.all(center[mask] >= 0)


def test_guest_graph_unwraps_a_molecule_longer_than_half_cell():
    atoms = [
        {"element": "C", "fractional": [x, 0.5, 0.5]}
        for x in [0.9, 0.05, 0.2, 0.35, 0.5]
    ]
    unwrapped = unwrap_group(list(range(5)), atoms, Lattice.cubic(10), {})
    assert np.allclose(np.diff(unwrapped[:, 0]), 0.15)
    assert unwrapped[-1, 0] - unwrapped[0, 0] > 0.5


def test_ordered_structure_required(tmp_path):
    path = tmp_path / "disordered.cif"
    Structure(Lattice.cubic(8), [{"C": 0.5}], [[0, 0, 0]]).to(filename=str(path))
    with pytest.raises(ValueError, match="partial occupancy"):
        intake(path)


def test_compute_display_and_cache_rejection(tmp_path):
    path = Path(__file__).resolve().parents[1] / "examples" / "demo-channel.cif"
    job = {"cif": str(path), "grid_step_A": 0.5, "pore_axis_uvw": [1, 0, 0]}
    calculation = compute(job, tmp_path / "calculation")
    geometry = display(job, tmp_path / "calculation", tmp_path / "display")
    assert calculation["passed"] and geometry["passed"]
    assert geometry["coordinate_modulo_lattice_error_A"] < 1e-10
    assert geometry["triangles"] > 100
    job["probe_radius_A"] = 1.3
    with pytest.raises(ValueError, match="cached"):
        display(job, tmp_path / "calculation", tmp_path / "other")


def test_invalid_job_and_resource_limit():
    with pytest.raises(ValueError, match="vdW"):
        radii({"radii_model": "custom", "vdw_radii_A": {"C": float("nan")}}, ["C"])
    with pytest.raises(ValueError, match="max_voxels"):
        fields(
            np.array([[0.0, 0.0, 0.0]]),
            ["C"],
            np.eye(3) * 20,
            {"C": 1.7},
            1.0,
            0.05,
            1000,
        )


def test_guest_removal_and_original_coordinate_display(tmp_path):
    path = Path(__file__).resolve().parents[1] / "examples" / "demo-channel.cif"
    source, _ = intake(path)
    source.append("C", [0.5, 0.5, 0.5])
    source.append("O", [0.5 - 1.16 / 12, 0.5, 0.5])
    source.append("O", [0.5 + 1.16 / 12, 0.5, 0.5])
    cif = tmp_path / "guest.cif"
    source.to(filename=str(cif))
    _, data = intake(cif)
    # Derive indices from this parsed structure, never from species ordering.
    ids = [
        a["source_index"]
        for a in data["atoms"]
        if np.allclose(a["fractional"][1:], [0.5, 0.5])
    ]
    ids.sort(key=lambda i: data["atoms"][i]["element"] != "C")
    job = {
        "cif": str(cif),
        "guest_groups": [ids],
        "grid_step_A": 0.5,
        "pore_axis_uvw": [1, 0, 0],
        "display_mode": "single",
        "crop_fractional": {"min": [0, 0.27, 0.27], "max": [1, 0.73, 0.73]},
    }
    result = compute(job, tmp_path / "calculation")
    geometry = display(job, tmp_path / "calculation", tmp_path / "display")
    assert result["framework_atoms"] == 128
    assert geometry["guest_molecules"] == 1
    output = json.loads((tmp_path / "display" / "display.json").read_text())
    assert len(output["guests"][0]["bonds"]) == 2
    assert geometry["coordinate_modulo_lattice_error_A"] < 1e-10
