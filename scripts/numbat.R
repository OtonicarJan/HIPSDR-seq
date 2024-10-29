library(data.table)
library(numbat)
library(Matrix)

# Prepare the reference
mtx <- read.csv("gene_mtx.csv", row.names=1)
mtx = t(mtx)
mat <- as(mtx, "dgCMatrix")
cells = read.delim("cell_types.tsv")
ref = aggregate_counts(mat, cells)

# Prepare sample, subset to matching RNA-DNA barcodes
mtx <- read.csv("hipsdr_mtx.csv", row.names=1)
mtx = t(mtx)
count_mat <- as(mtx, "dgCMatrix")
df_allele <- fread("lfs_hipsdr_allele_counts.tsv.gz")

out <- run_numbat(
    count_mat=count_mat,
    lambdas_ref=ref,
    df_allele=df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = 8,
    plot = TRUE,
    out_dir = './HIPSDR_fibroblasts'
)