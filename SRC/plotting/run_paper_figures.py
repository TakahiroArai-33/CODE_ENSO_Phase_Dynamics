#!/usr/bin/env python3
"""Generate the selected ERSSTv5 paper panels and publish complete PNG/EPS pairs."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = REPO_ROOT / "SRC" / "enso_ERSSTv5_2026"
OUTPUT_DIR = REPO_ROOT / "SRC" / "figure"
PLOT_SCRIPTS = ("run_selected_panels.py",)
MANUSCRIPT_REFERENCE = "Arai and Kawamura, ENSO phase dynamics manuscript (2026-10-30)"


def panel_mapping() -> list[tuple[str, str]]:
    """The 44 constituent panels used in the current manuscript, in paper order."""
    panels: list[tuple[str, str]] = []

    def add(source: str, folder: str, suffix: str) -> None:
        panels.append((source, f"{folder}/{folder}_{suffix}"))

    for mode, name, letter in ((2, "Annual", "a"), (9, "ENSO", "b")):
        add(f"ComplexPlane_g{mode}_amp30_fixed", "fig1", f"{letter}_complex_plane_{name.lower()}")
        for phase in range(1, 9):
            add(f"PhaseComposite_{name}_amp30_fixed_j{phase}", "fig1",
                f"{letter}_phase_composite_{name.lower()}_j{phase}")

    add("Q_e_amp30_M3_ERSSTv5_2026", "fig2", "a_Qe")
    add("Q_e_ODE_amp30_M3_ERSSTv5_2026_ODEoverlay", "fig2", "b_Qe_ODE_overlay")
    for event, letter in (("ElNino", "c"), ("LaNina", "d")):
        add(f"Q_e_passage_{event}_amp30_M3_ERSSTv5_2026_SDE", "fig2",
            f"{letter}_{event.lower()}_duration_SDE")

    avg = "Averaging_amp30_M3_ratio1to3_ERSSTv5_2026_psi_cycles"
    direct = "DirectAveraging_amp30_M3_ratio1to3_ERSSTv5_2026_psi_cycles"
    add(f"{avg}_Psi", "fig3", "psi_ratio3to1")
    add(f"{direct}_Psi", "fig4", "a_psi_ratio3to1")
    add(f"{direct}_Hist", "fig4", "b_hist_ratio3to1")
    add("nino_reconstruction_double_peaks", "fig5", "nino34_reconstruction_double_peaks")

    for event, month, letter in (("ElNino", 4, "a"), ("ElNino", 9, "b"),
                                 ("LaNina", 4, "c"), ("LaNina", 9, "d")):
        add(f"Langevin_passage_{event}_amp30_M3_ERSSTv5_2026_estimated_month{month:02d}",
            "fig6", f"{letter}_{event.lower()}_month{month:02d}")
    for letter, amp in zip("abcde", (10, 20, 30, 40, 50)):
        add(f"Q_e_amp{amp}_M3_ERSSTv5_2026_threshold_title", "fig_s1", f"{letter}_Qe_amp{amp}")
    add("DirectAveraging_overlay_M3_ratio1to3_ERSSTv5_2026_psi_cycles_Psi",
        "fig_s2", "psi_ratio3to1_amp10to50")
    for letter, source_ratio, paper_ratio in (("a", "1to2", "2to1"), ("b", "1to4", "4to1"),
                                             ("c", "1to5", "5to1"), ("d", "2to5", "5to2")):
        for panel in ("Psi", "Hist"):
            add(f"DirectAveraging_amp30_M3_ratio{source_ratio}_ERSSTv5_2026_psi_cycles_{panel}",
                "fig_s3", f"{letter}_{panel.lower()}_ratio{paper_ratio}")
    return panels


def required_inputs(include_estimates: bool = True) -> list[Path]:
    src = DATA_DIR / "SRC_py"
    paths = [DATA_DIR / "myFolder" / name for name in
             ("Z.mat", "Z_left.mat", "mu.mat", "Lambda.mat", "Nino34_data.mat")]
    paths.append(DATA_DIR / "CLIMATE_DATA" / "sst.mnmean.nc")
    if not include_estimates:
        return paths
    paths.extend(src / "Interim_data_Q" / f"Q_result_fixed_M3_amp{amp}.npz"
                 for amp in (10, 20, 30, 40, 50))
    paths.append(src / "Interim_Averaging" / "Q_result_fixed_M3_amp30_ratio1:3_averaged_result.npz")
    direct_pairs = [(30, ratio) for ratio in ("1:2", "1:3", "1:4", "1:5", "2:5")]
    direct_pairs.extend((amp, "1:3") for amp in (10, 20, 40, 50))
    paths.extend(src / "Interim_data_DirectAveraging" / f"DirectAvaraging_result_M3_amp{amp}_ratio{ratio}.npz"
                 for amp, ratio in direct_pairs)
    return paths


def collect_sources(source_dir: Path) -> list[tuple[Path, str]]:
    """Validate requested pairs, including outputs indexed by an existing manifest."""
    source_dir = source_dir.resolve()
    wanted = {f"{stem}.{extension}" for stem, _ in panel_mapping() for extension in ("png", "eps")}
    found: dict[str, Path] = {}

    def register(name: str, path: Path) -> None:
        if name not in wanted or not path.is_file():
            return
        if not path.stat().st_size:
            raise ValueError(f"Empty figure: {path}")
        if name in found and found[name] != path:
            if hashlib.sha256(found[name].read_bytes()).digest() != hashlib.sha256(path.read_bytes()).digest():
                raise ValueError(f"Conflicting generated filename: {found[name]} and {path}")
        else:
            found[name] = path

    for path in source_dir.rglob("*"):
        register(path.name, path)
    manifest = source_dir / "figure_manifest.tsv"
    if manifest.is_file():
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = (source_dir / row["file"]).resolve()
                if not path.is_relative_to(source_dir):
                    raise ValueError(f"Manifest file is outside the source directory: {path}")
                register(row["generated_basename"], path)
    missing = sorted(wanted - found.keys())
    if missing:
        raise FileNotFoundError("Missing paper figures:\n" + "\n".join(missing))
    return [(found[f"{stem}.{extension}"], f"{target}.{extension}")
            for stem, target in panel_mapping() for extension in ("png", "eps")]


def publish(source_dir: Path, output_dir: Path, metadata: dict, dry_run: bool = False) -> int:
    sources = collect_sources(source_dir)
    if dry_run:
        for source, target in sources:
            print(f"{source.name} -> {target}")
        return len(sources)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".paper-publish-", dir=output_dir.parent) as work:
        work_dir = Path(work)
        selected = work_dir / "selected"
        selected.mkdir()
        manifest = ["file\tgenerated_basename"]
        generated_names = {
            f"{target}.{extension}": f"{stem}.{extension}"
            for stem, target in panel_mapping() for extension in ("png", "eps")
        }
        for source, target in sources:
            destination = selected / target
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            manifest.append(f"{target}\t{generated_names[target]}")
        (selected / "figure_manifest.tsv").write_text("\n".join(manifest) + "\n")
        metadata = {**metadata, "figure_files": len(sources)}
        (selected / "run_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
        previous = work_dir / "previous"
        if output_dir.exists():
            output_dir.rename(previous)
        try:
            selected.rename(output_dir)
        except BaseException:
            if previous.exists():
                previous.rename(output_dir)
            raise
    return len(sources)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--collect-from", type=Path, help="Collect previously generated PNG/EPS pairs without plotting.")
    parser.add_argument("--dry-run", action="store_true", help="Validate a --collect-from directory without copying.")
    parser.add_argument("--preflight-estimation", action="store_true", help="Check inputs and options before estimation, then exit.")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir == SCRIPT_DIR or output_dir in SCRIPT_DIR.parents:
        parser.error("output-dir must be a dedicated generated-output directory")
    if (output_dir.exists() and output_dir != OUTPUT_DIR and any(output_dir.iterdir())
            and not (output_dir / "run_metadata.json").is_file()):
        parser.error("A nonempty custom output-dir must contain a previous run_metadata.json")
    if args.dry_run and args.collect_from is None:
        parser.error("--dry-run requires --collect-from")
    if args.preflight_estimation and (args.collect_from is not None or args.dry_run):
        parser.error("collection options cannot be used with --preflight-estimation")
    if args.collect_from is not None:
        source_dir = args.collect_from.resolve()
        previous_metadata = source_dir / "run_metadata.json"
        if previous_metadata.is_file():
            metadata = json.loads(previous_metadata.read_text())
            metadata["previous_mode"] = metadata.get("mode")
            metadata["mode"] = "reselect_existing"
        else:
            metadata = {"mode": "collect_existing", "sde_trial_count": None,
                        "note": "Simulation settings of existing figures were not verified."}
        metadata["selection_reference"] = MANUSCRIPT_REFERENCE
        metadata["selection_updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        count = publish(source_dir, output_dir, metadata, args.dry_run)
        print(f"Paper figure files: {count} (dry-run={args.dry_run})")
        return 0

    paths = required_inputs(include_estimates=not args.preflight_estimation)
    missing = [str(p.relative_to(REPO_ROOT)) for p in paths if not p.is_file() or not p.stat().st_size]
    if missing:
        parser.error("Missing inputs. Download NOAA SST with bash SRC/download_sst_data.sh; "
                     "create estimates with the dataset run_all.sh:\n" + "\n".join(missing))
    missing_modules = [name for name in ("numpy", "scipy", "matplotlib", "netCDF4", "pandas")
                       if importlib.util.find_spec(name) is None]
    if missing_modules:
        parser.error("Missing Python packages: " + ", ".join(missing_modules))
    missing_scripts = [name for name in PLOT_SCRIPTS if not (SCRIPT_DIR / name).is_file()]
    if missing_scripts:
        parser.error("Missing plotting scripts: " + ", ".join(missing_scripts))
    try:
        trials = int(os.environ.get("SDE_TRIAL_COUNT", "5000"))
    except ValueError:
        parser.error("SDE_TRIAL_COUNT must be a positive integer")
    if trials < 1:
        parser.error("SDE_TRIAL_COUNT must be a positive integer")
    if args.preflight_estimation:
        print("Input, environment, and figure options checked.")
        return 0
    metadata = {
        "selection_reference": MANUSCRIPT_REFERENCE,
        "mode": "generate", "dataset": "ERSSTv5_2026", "Q_M": 3,
        "sde_trial_count": trials, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "scripts": list(PLOT_SCRIPTS),
    }
    with tempfile.TemporaryDirectory(prefix=".paper-figure-stage-", dir=SCRIPT_DIR) as work:
        env = {**os.environ, "PAPER_FIGURE_OUTPUT_ROOT": work, "MPLBACKEND": "Agg",
               "SDE_TRIAL_COUNT": str(trials)}
        for index, name in enumerate(PLOT_SCRIPTS, 1):
            started = time.monotonic()
            print(f"[{index}/{len(PLOT_SCRIPTS)}] {name}", flush=True)
            subprocess.run([sys.executable, str(SCRIPT_DIR / name)], cwd=REPO_ROOT, env=env, check=True)
            print(f"  finished in {time.monotonic() - started:.1f} s", flush=True)
        count = publish(Path(work), output_dir, metadata)
    print(f"Saved {count} paper figure files to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
