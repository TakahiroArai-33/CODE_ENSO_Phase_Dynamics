#!/usr/bin/env python3
# coding: utf-8
"""Render only panels used in ENSO_Fig1–6 and ENSO_FigS1–S3.

The publication runner sets PAPER_FIGURE_OUTPUT_ROOT to its staging directory
and collects these 44 PNG/EPS pairs into the numbered figure directories.
Each plotting script runs in a separate process to retain its plotting style.
"""

from pathlib import Path
import subprocess
import sys


PLOT_SCRIPTS = (
    "PhaseComposite.py",                       # Fig1: 18 panels
    "Plot_Qe.py",                              # Fig2: 1; FigS1: 5
    "Plot_Qe_ODEOverlay.py",                   # Fig2: 1
    "Q_eの通過時間計算.py",                      # Fig2: 1
    "Q_eの通過時間計算_ラニーニャ.py",            # Fig2: 1
    "Plot_Averaging_PsiCycles.py",              # Fig3: 1
    "Plot_DirectAveraging_PsiCycles.py",         # Fig4: 2; FigS3: 8
    "NinoSeries.py",                           # Fig5: 1
    "Plot_Langevin_Passage_Cartoon.py",          # Fig6: 4
    "Plot_DirectAveraging_AmpOverlay_PsiCycles.py",  # FigS2: 1
)


def main():
    script_dir = Path(__file__).resolve().parent
    for script_name in PLOT_SCRIPTS:
        print(f"Rendering {script_name}", flush=True)
        subprocess.run([sys.executable, str(script_dir / script_name)],
                       cwd=script_dir, check=True)


if __name__ == "__main__":
    main()
