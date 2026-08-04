"""CLI argument parsing tests for sar-sea-ice-mapping."""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from conftest import SCRIPT


def run_cli(args, timeout=90):
    return subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True, timeout=timeout,
    )


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--season" in r.stdout
        assert "--synthetic" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        r = run_cli([])
        assert r.returncode == 2

    def test_bad_season_exit_2(self):
        r = run_cli(["--bbox", "120", "75", "122", "77", "--season", "spring"])
        assert r.returncode == 2

    def test_input_not_found_exit_2(self):
        r = run_cli(["--input", "nonexistent.tif"])
        assert r.returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        r = run_cli(["--synthetic"])
        assert r.returncode == 2


class TestCLIEndToEnd:
    def test_bbox_only_runs(self, tmp_path):
        out = str(tmp_path / "o1")
        r = run_cli(["--bbox", "120", "75", "122", "77", "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "ice_type.tif"))
        assert os.path.exists(os.path.join(out, "ice_concentration.tif"))
        assert os.path.exists(os.path.join(out, "ice_statistics.json"))
        assert os.path.exists(os.path.join(out, "output-manifest.json"))

    def test_summer_season(self, tmp_path):
        out = str(tmp_path / "o2")
        r = run_cli([
            "--bbox", "120", "75", "122", "77", "--synthetic",
            "--season", "summer", "--output-dir", out, "--quiet",
        ])
        assert r.returncode == 0, r.stderr

    def test_detects_ice_types(self, tmp_path):
        """合成场景应分出冰类（含多年冰）。"""
        import json
        out = str(tmp_path / "o3")
        r = run_cli(["--bbox", "120", "75", "122", "77", "--synthetic",
                     "--output-dir", out, "--quiet"])
        assert r.returncode == 0, r.stderr
        with open(os.path.join(out, "ice_statistics.json"), encoding="utf-8") as f:
            stats = json.load(f)
        assert stats["ice_fraction"] > 0.1
        assert stats["per_class"]["multiyear_ice"]["pixels"] > 0
