#!/usr/bin/env Rscript
# Export a Seurat .rds to a plain on-disk bundle that Python can assemble
# into an h5ad, without SeuratDisk.
#
# 06_2 used SaveH5Seurat() + Convert(). SeuratDisk is unmaintained and its
# Convert() frequently fails on Seurat v5 objects (the v5 Assay class stores
# counts in layers, which the h5Seurat writer does not understand). Writing
# Matrix Market + CSV avoids the dependency entirely: the formats are stable,
# the failure modes are obvious, and the result is assembled in Python by
# assemble_h5ad.py.
#
# Usage:
#   Rscript convert_seurat_to_h5ad.R <input.rds> <outdir> [assay]
#
# Writes into <outdir>/:
#   counts.mtx        genes x cells raw counts (Matrix Market)
#   barcodes.csv      cell names, in matrix column order
#   features.csv      gene names, in matrix row order
#   metadata.csv      the full meta.data table
#   reduction_<k>.csv one file per dimensionality reduction
#   manifest.json     shapes, assay name, layer used, reduction list

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript convert_seurat_to_h5ad.R <input.rds> <outdir> [assay]")
}
in_rds <- args[[1]]
outdir <- args[[2]]
assay_arg <- if (length(args) >= 3) args[[3]] else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) {
  stop(sprintf("%s contains a %s, not a Seurat object", in_rds, class(obj)[1]))
}

assay <- if (is.null(assay_arg)) DefaultAssay(obj) else assay_arg
message(sprintf("object: %d features x %d cells | assays: %s | using: %s",
                nrow(obj), ncol(obj), paste(Assays(obj), collapse = ","), assay))

# Prefer raw counts; fall back to data when an object ships normalised only.
layer_used <- "counts"
mat <- tryCatch(
  SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
  error = function(e) NULL
)
if (is.null(mat) || length(mat) == 0 || nrow(mat) == 0) {
  message("no counts layer; falling back to 'data' (values will not be integers)")
  layer_used <- "data"
  mat <- SeuratObject::LayerData(obj, assay = assay, layer = "data")
}
mat <- as(mat, "CsparseMatrix")

invisible(writeMM(mat, file.path(outdir, "counts.mtx")))
write.csv(data.frame(barcode = colnames(mat)),
          file.path(outdir, "barcodes.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(feature = rownames(mat)),
          file.path(outdir, "features.csv"), row.names = FALSE, quote = TRUE)

meta <- obj@meta.data
meta <- meta[colnames(mat), , drop = FALSE]   # enforce matrix column order
write.csv(meta, file.path(outdir, "metadata.csv"), row.names = TRUE, quote = TRUE)

reductions <- Reductions(obj)
for (r in reductions) {
  emb <- Embeddings(obj, reduction = r)
  emb <- emb[colnames(mat), , drop = FALSE]
  write.csv(emb, file.path(outdir, sprintf("reduction_%s.csv", r)),
            row.names = TRUE, quote = TRUE)
}
message(sprintf("reductions written: %s",
                if (length(reductions)) paste(reductions, collapse = ", ") else "none"))

manifest <- sprintf(
  '{"source":"%s","assay":"%s","layer":"%s","n_features":%d,"n_cells":%d,"reductions":[%s],"meta_columns":%d}',
  basename(in_rds), assay, layer_used, nrow(mat), ncol(mat),
  if (length(reductions)) paste(sprintf('"%s"', reductions), collapse = ",") else "",
  ncol(meta)
)
writeLines(manifest, file.path(outdir, "manifest.json"))
message(sprintf("wrote bundle to %s", outdir))
