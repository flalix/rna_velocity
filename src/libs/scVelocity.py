"""Load kb-python nac output for the three K409 libraries, call cells, save h5ad.

Layers: spliced = mature, unspliced = nascent, ambiguous kept separately.
X = mature + ambiguous (cytoplasmic convention, closest to CellRanger).

Expects, relative to the working directory (override the root with self.COUNT_ROOT):
    count_<library>/counts_unfiltered/cells_x_genes.{mature,nascent,ambiguous}.mtx
    geo_barcodes/author_cells.json      <- run fetch_author_barcodes.py first
"""
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

class Velocitye(object):
    def __init__(self, tumor_primary:str='SRR11492025', lymph_node:str='SRR11492017', blood_PBMC:str='SRR11492013'):

        self.COUNT_ROOT = os.environ.get("self.COUNT_ROOT", ".")
        self.BARCODE_JSON = os.environ.get("AUTHOR_CELLS", "geo_barcodes/author_cells.json")

        # tumor_primary:SRR11492025 lymph_node:SRR11492017 blood_PBMC:SRR11492013
        self.LIBS = {"tumor_primary": tumor_primary,
                     "lymph_node": "involved regional lymph node",
                     "blood_PBMC": "peripheral blood"}

    def load_one(self, lib):
        d = os.path.join(self.COUNT_ROOT, f"count_{lib}", "counts_unfiltered")
        bc = pd.read_csv(f"{d}/cells_x_genes.barcodes.txt", header=None)[0].astype(str).values
        genes = pd.read_csv(f"{d}/cells_x_genes.genes.txt", header=None)[0].astype(str).values
        names_f = f"{d}/cells_x_genes.genes.names.txt"
        names = pd.read_csv(names_f, header=None)[0].astype(str).values if os.path.exists(names_f) else genes
        layers = {k: sc.read_mtx(f"{d}/cells_x_genes.{k}.mtx").X.astype(np.float32)
                for k in ("mature", "nascent", "ambiguous")}
        A = ad.AnnData(X=(layers["mature"] + layers["ambiguous"]).tocsr(),
                    obs=pd.DataFrame(index=[f"{lib}:{b}" for b in bc]),
                    var=pd.DataFrame({"gene_id": genes, "gene_name": names}, index=names))
        A.layers["spliced"] = layers["mature"].tocsr()
        A.layers["unspliced"] = layers["nascent"].tocsr()
        A.layers["ambiguous"] = layers["ambiguous"].tocsr()
        A.obs["library"] = lib
        A.obs["tissue"] = self.LIBS[lib]
        A.obs["barcode"] = bc
        A.var_names_make_unique()
        return A


    def knee_threshold(counts, floor=500):
        """Inflection of the log-log barcode-rank curve (kneedle-style), with a UMI floor."""
        v = np.sort(counts)[::-1]
        v = v[v > 0]
        x, y = np.log10(np.arange(1, len(v) + 1)), np.log10(v)
        x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
        num = np.abs((y1 - y0) * x - (x1 - x0) * y + (x1 * y0 - y1 * x0))
        idx = int(np.argmax(num / np.hypot(y1 - y0, x1 - x0)))
        return max(float(v[idx]), floor), idx, v

