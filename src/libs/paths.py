"""Project paths, resolved from this file's location.

Avoids hard-coded absolute paths and avoids the Claude Science artifact
store, whose layout (.../artifacts/<project>/<artifact_id>/v<prefix>_<name>)
is internal and not a stable location to reference from analysis code.

    from libs.paths import PROJ, H5AD, TABLES, FIGURES
    IN = sc.read_h5ad(H5AD / "innervation_selection_scfates.h5ad")
    mods = pd.read_csv(TABLES / "innervation_gene_modules.csv", index_col=0)
"""

from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]      # src/libs/paths.py -> project root
DATA = PROJ / "data"
H5AD = DATA / "h5ad"
BUNDLES = DATA / "bundles"
RESULTS = PROJ / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

for _d in (H5AD, BUNDLES, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
