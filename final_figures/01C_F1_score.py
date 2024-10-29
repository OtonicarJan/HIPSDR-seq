import argparse
import logging

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score


def arguments():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description="Bin size")
    parser.add_argument("bin", type=int, help="Bin size")
    args = parser.parse_args()

    return args


args = arguments()
binsize = args.bin
# Set up basic configuration for logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Create logger object
logger = logging.getLogger(__name__)

BIN_SIZE = [binsize]

cl0_rna = []
with open("../data/scrna_cluster0.txt") as fin:
    for line in fin:
        cl0_rna.append(line.rstrip())

cl0_dna = []
with open("../data/HIPSDR_cluster0.txt") as fin:
    for line in fin:
        cl0_dna.append(line.rstrip())

metacells_cl0 = []
with open("../data/metacells_cluster0.txt", "r") as fin:
    for line in fin:
        metacells_cl0.append(line.rstrip())


def compute_f1(bulk, single_cell, comparison):
    res = []
    for bs in [5, 10, 50, 100]:
        logger.info(f"Working on {bs}")
        for i in range(0, len(bulk), bs):
            logger.info(f"Working on {i}")
            start = i
            end = i + bs
            group = single_cell.iloc[:, start:end]
            short_bulk = bulk[start:end]
            # only analyse window if not all bins are equal to 2
            if not (short_bulk == 2).all():
                logger.info(f"Computing F1 scores")
                f1_scr = group.apply(
                    lambda row: f1_score(row, short_bulk, average="macro"), axis=1
                )
                avg_score = f1_scr.mean()
                res.append((avg_score, comparison, bs, i))

    res = pd.DataFrame(
        res, columns=["f1_score", "relation", "window_size", "bin"]
    )
    res["window_size"] = pd.Categorical(
        res["window_size"], categories=[5, 10, 50, 100], ordered=True
    )

    return res


def prepare_bulk_for_dna(path):
    bulk = pd.read_csv(path, sep="\t", header=None)
    bulk.columns = ["chrom_num", "start", "end", "CN", "median"]
    bulk.chrom_num = [x[3:] for x in bulk.chrom_num]

    bulk["chrom"] = "chr" + bulk.chrom_num
    bulk.chrom_num = bulk.chrom_num.replace({"X": 23, "Y": 24})
    bulk.chrom_num = bulk.chrom_num.astype(int)
    bulk = bulk.sort_values(by=["chrom_num", "start"])
    bulk["CN"] = bulk["CN"] - 1
    bulk["CN"] = bulk["CN"].replace(0, 1)
    bulk["CN"] = np.where(bulk["CN"] < 2, 1, np.where(bulk["CN"] > 2, 3, bulk["CN"]))

    bulk["window"] = [f"{tup.chrom}:{tup.start}-{tup.end}" for tup in bulk.itertuples()]
    bulk = bulk.loc[bulk["chrom"] != "chrY"]
    bulk.set_index("window", inplace=True)
    return bulk


def categorize_row(row, first_cutoff, second_cutoff):
    # Calculate the cutoffs for loss/gain
    loss_cutoff = np.quantile(row, first_cutoff)
    gain_cutoff = np.quantile(row, second_cutoff)

    # Apply the transformation based on the terciles
    return np.where(row < loss_cutoff, 1, np.where(row > gain_cutoff, 3, 2))


def create_custom_dataframe(n, num_columns, x_percent, y_percent, z_percent):
    data = []

    # Calculate the number of 1s, 2s, and 3s based on the percentages
    num_ones = int(num_columns * (x_percent / 100))
    num_twos = int(num_columns * (y_percent / 100))
    num_threes = num_columns - (num_ones + num_twos)

    for _ in range(n):
        # Create a row with the correct number of 1s, 2s, and 3s
        row = [1] * num_ones + [2] * num_twos + [3] * num_threes
        np.random.shuffle(row)  # Shuffle the row to randomize the order

        data.append(row)

    df = pd.DataFrame(data)

    return df


for bin_size in BIN_SIZE:
    logger.info(f"Working on bin {bin_size}")
    p195 = prepare_bulk_for_dna(
        f"../data/readcounts.{bin_size}kb.cell_bc_p195_per_bin_calls_long.bed"
    )

    single_dna = pd.read_csv(
        f"../data/CNVs_HIPSDR_filtered_{bin_size}kb.csv", index_col=0
    )
    single_dna = single_dna.loc[[cl for cl in cl0_dna if cl in single_dna.index]]

    metacells_dna = pd.read_csv(f"../data/cnvs_metacells_{bin_size}kb.csv", index_col=0)
    metacells_dna = metacells_dna.loc[metacells_cl0]

    single_rna = pd.read_csv(
        f"../data/numbat_rna_windows_{bin_size}.csv.gz", index_col=0
    )
    single_rna = single_rna.loc[cl0_rna]

    # Categorize smoothed expression values to either loss, neutral or gain (1, 2, 3)
    first_tercile = round(
        p195.loc[single_rna.columns]["CN"].value_counts()[1]
        / len(p195.loc[single_rna.columns]),
        3,
    )
    second_tercile = 1 - round(
        p195.loc[single_rna.columns]["CN"].value_counts()[3]
        / len(p195.loc[single_rna.columns]),
        3,
    )
    single_rna_categorized = single_rna.apply(
        categorize_row, axis=1, args=(first_tercile, second_tercile)
    )

    single_rna = pd.DataFrame(
        single_rna_categorized.to_list(),
        columns=single_rna.columns,
        index=single_rna.index,
    )

    overlap_cols = single_rna.columns

    single_dna = single_dna.loc[:, overlap_cols].copy()
    single_dna = single_dna.map(lambda x: 3 if x > 3 else x)
    metacells_dna = metacells_dna.loc[:, overlap_cols].copy()
    metacells_dna = metacells_dna.map(lambda x: 3 if x > 3 else x)
    random_control = create_custom_dataframe(
        500,
        single_dna.shape[1],
        p195.loc[single_rna.columns]["CN"].value_counts()[1]
        / len(p195.loc[single_rna.columns])
        * 100,
        p195.loc[single_rna.columns]["CN"].value_counts()[2]
        / len(p195.loc[single_rna.columns])
        * 100,
        p195.loc[single_rna.columns]["CN"].value_counts()[3]
        / len(p195.loc[single_rna.columns])
        * 100,
    )
    diploid_control = pd.DataFrame(np.full((2, single_dna.shape[1]), 2))
    bulk = p195.loc[overlap_cols]["CN"]

    f_ones_dna = compute_f1(bulk, single_dna, "bulk-scDNA")
    f_ones_metacells = compute_f1(bulk, metacells_dna, "bulk-metacells")
    f_ones_rna = compute_f1(bulk, single_rna, "bulk-scRNA")
    f_ones_control = compute_f1(bulk, random_control, "bulk-random_control")
    f_ones_diploid = compute_f1(bulk, diploid_control, "bulk-diploid_control")

    merged = pd.concat(
        [f_ones_metacells, f_ones_dna, f_ones_rna, f_ones_control, f_ones_diploid]
    )

    merged.to_csv(f"merged_bin_{bin_size}kb.csv")
