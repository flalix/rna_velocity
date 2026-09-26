"""First cell for any notebook in this project.

    exec(open("_header.py").read())        # or paste the body into cell 1

Does the three things a notebook actually needs, none of which os.system can
do (a child shell cannot change the parent kernel's cwd, environment or
interpreter):

1. puts src/ on sys.path, found by walking up from the notebook's directory,
   so the same cell works from notebooks/, the project root, or a subfolder
2. imports libs.run_scfates FIRST, which resolves R_HOME and puts R on PATH
   before scFates is imported -- both are cached at import and cannot be
   fixed afterwards
3. prints the three facts worth checking before any analysis

The interpreter comes from the kernel picker, not from code: select the
.venv kernel. `conda activate` and `.venv/bin/activate` have no effect from
inside a running kernel.
"""

import os
import shutil
import sys
from pathlib import Path

SRC = next(p / "src" for p in [Path.cwd(), *Path.cwd().parents]
           if (p / "src" / "libs").is_dir())
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libs.paths import PROJ, DATA, H5AD, BUNDLES, RESULTS, TABLES, FIGURES
from libs.run_scfates import preprocess, fit_curve, plot_trajectory, test_association

print("python :", sys.executable)
print("PROJ   :", PROJ)
print("R      :", shutil.which("R"), "| R_HOME:", os.environ.get("R_HOME"))
