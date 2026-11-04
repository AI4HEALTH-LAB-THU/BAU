"""Input and output paths for the analyses."""

import os
from pathlib import Path

CODE_ROOT = Path(os.environ.get("BAU_CODE_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = Path(os.environ.get("BAU_DATA_ROOT", "/home/dnanexus/bau-work/data")).resolve()
RESULTS_ROOT = Path(os.environ.get("BAU_RESULTS_ROOT", "/home/dnanexus/bau-work/results")).resolve()
