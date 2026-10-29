"""Portable cache and figure destinations for the paper plotting scripts."""
from pathlib import Path
from tempfile import gettempdir
import os

_CACHE_DIR = Path(gettempdir()) / "ersstv5-paper-cache"
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_DIR))
os.environ.setdefault("MPLBACKEND", "Agg")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
FIGURE_OUTPUT_DIR = Path(os.environ.get(
    "PAPER_FIGURE_OUTPUT_ROOT", str(Path(__file__).resolve().parent.parent / "figure" / "materials")
)).resolve()
FIGURE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
