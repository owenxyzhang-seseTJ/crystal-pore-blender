from pathlib import Path
import subprocess
import sys

from crystal_pore_blender.cli import main


def test_module_propagates_failure_status():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "crystal_pore_blender",
            "run",
            "--out",
            "unused",
            "--no-render",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "exactly one" in result.stderr


def test_init_and_validate(tmp_path, monkeypatch):
    cif = Path(__file__).resolve().parents[1] / "examples" / "demo-channel.cif"
    target = tmp_path / "job.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["crystal-pore-blender", "init", "--cif", str(cif), "--out", str(target)],
    )
    assert main() == 0
    monkeypatch.setattr(
        sys, "argv", ["crystal-pore-blender", "validate", "--job", str(target)]
    )
    assert main() == 0
    assert target.with_name("job.intake.json").exists()
