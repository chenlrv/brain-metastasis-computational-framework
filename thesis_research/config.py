"""Location of the project data, read from environment variables.

THESIS_PROJECT_ROOT (required)
    Folder holding the data and results, laid out as described in the README
    (resources/, outputs/cell_annotation/, thesis_plots/, automation/outputs/).
    Every analysis script reads its inputs from, and writes its outputs to,
    paths under this folder. There is deliberately no default, so that no
    script can write anywhere the user has not chosen.

THESIS_L321_TX_FILE (FastReseg steps only)
    The raw transcript table of slide L321 (L321_tx_file.csv) from the
    vendor's flat-file export.
"""
import os
import pathlib

_root = os.environ.get("THESIS_PROJECT_ROOT", "").strip()
if not _root:
    raise RuntimeError(
        "THESIS_PROJECT_ROOT is not set. Set it to the folder holding the project "
        "data (see README, section 'Data layout') before running any script.")

PROJECT_ROOT = pathlib.Path(_root)
PROJECT_ROOT_STR = PROJECT_ROOT.as_posix()
L321_TX_FILE = os.environ.get("THESIS_L321_TX_FILE", "")
