"""Transcript reassignment test, version 3: parameters tied to cell size, no current-cell advantage.

Runs reassign_v2.py's procedure unchanged, with parameters that each have a
stated reason instead of being round numbers in pixel or log units:
  * candidate radius 14.4 um   -- slightly more than one median cell diameter (12.6 um)
  * up to 5 candidate cells    -- allowing up to 20 changed no result (reassign_v2_k_sweep.py)
  * decay length 6.3 um        -- the median cell radius, sqrt(123.8 um^2 / pi)
  * no current-cell advantage  -- the test is meant to be permissive: each transcript
                                  goes to whichever candidate explains it best
Sweeps are centred on these values. Writes only to
agents/outputs/segmentation_reassign_v3/ (must not exist); earlier runs are untouched.

Run: conda run -n thesis_research python agents/segmentation/reassign_v3.py
"""
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import reassign_v2 as rv  # noqa: E402
import reassign_v2_k5 as k5  # noqa: E402  (provides own_cell_coverage; sets K_NEIGH = 5)

CELL_RADIUS_UM = float(np.sqrt(123.8 / np.pi))          # 6.28 um, median cell area in slice 1
rv.K_NEIGH = 5
rv.OUT = rv.ROOT / "agents/outputs/segmentation_reassign_v3"
rv.CANON = {"radius_um": 14.4, "decay_um": round(CELL_RADIUS_UM, 1), "home_adv": 0.0}
rv.SWEEPS = {
    "radius_um": [7.2, 14.4, 24.1],
    "home_adv": [np.inf, 3.0, 1.5, 0.5, 0.0],
    "decay_um": [round(CELL_RADIUS_UM * f, 1) for f in (0.25, 0.5, 1, 2)] + [np.inf],
}

if __name__ == "__main__":
    print("canonical:", rv.CANON, "| decay sweep:", rv.SWEEPS["decay_um"])
    rv.main()                                            # refuses if the output folder exists
    cov = k5.own_cell_coverage()
    json.dump(cov, open(rv.OUT / "own_cell_coverage.json", "w"), indent=2)
    print(f"own cell among candidates: {cov['pct_own_cell_candidate']:.2f}%")
