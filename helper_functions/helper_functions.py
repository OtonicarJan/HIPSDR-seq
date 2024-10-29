import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from natsort import natsorted
from tqdm import tqdm


def prepare_counts_df(path, metacelling=False, binsize=100000):
    barcodes = []
    step = binsize
    counts = []
    for i, file_name in enumerate(tqdm(path)):
        cs = []
        try:
            file = pd.read_csv(file_name, sep="\t", header=None)
            barcode = file_name.split("/")[-1].split("_")[-1].split(".")[0]
            barcodes.append(barcode)
            for tup in file.itertuples():
                try:
                    c = int(tup[1])
                    cs.append(c)
                except ValueError:
                    pass
            counts.append(cs)
        except:
            print(file_name)
    for i, file_name in enumerate(tqdm(path)):
        chrom, start = ["", ""]

        if i == 0:
            regions = []
            print("recording")
            file = pd.read_csv(file_name, sep="\t", header=None)
            for tup in file.itertuples():
                try:
                    c = int(tup[1])
                    regions.append((chrom, start0))
                    start0 = start0 + step
                except ValueError:
                    chrom, start, _ = [x.split("=")[1] for x in tup[1].split(" ")[1:-1]]
                    start0 = int(start)
                    pass

    counts_df_orig = pd.DataFrame(counts)
    counts_df_orig.index = barcodes
    counts_df_orig.columns = [f"{x[0]}:{x[1]}-{x[1]+step}" for x in regions]

    if metacelling:
        return counts_df_orig, regions
    else:
        return counts_df_orig


def prepare_cnvs(bed_files):
    cnas = []

    for file in tqdm(bed_files):
        if file.endswith(".bed"):
            try:
                cell_file = pd.read_csv(file, header=None, sep="\t")

                cell_file["bin"] = (
                    cell_file[0].astype(str)
                    + ":"
                    + cell_file[1].astype(str)
                    + "-"
                    + cell_file[2].astype(str)
                )
                cell_file = cell_file.set_index("bin")
                cell_file = cell_file[[3]]
                cnas.append(cell_file)
                cell = file.split("/")[-1].split("_")[2]
                cell_file.columns = [cell]
            except:
                print(file)
                continue

    cna = pd.concat(cnas, ignore_index=False, axis=1)
    cna = cna.T
    cna = cna.reindex(columns=natsorted(cna.columns))

    return cna


def lorenz(arr):
    # this divides the prefix sum by the total sum
    # this ensures all the values are between 0 and 1.0
    arr = np.sort(arr)
    scaled_prefix_sum = arr.cumsum() / arr.sum()
    # this prepends the 0 value (because 0% of all people have 0% of all wealth)
    return np.insert(scaled_prefix_sum, 0, 0)


def gini(arr):
    sorted_arr = arr.copy()
    sorted_arr = np.sort(sorted_arr)
    n = arr.size
    coef_ = 2.0 / n
    const_ = (n + 1.0) / n
    weighted_sum = sum([(i + 1) * yi for i, yi in enumerate(sorted_arr)])
    return coef_ * weighted_sum / (sorted_arr.sum()) - const_


def prepare_bulk(path):
    bulk = pd.read_csv(path, sep="\t", header=None)
    bulk.columns = ["chrom_num", "start", "end", "CN", "median"]
    bulk.chrom_num = [x[3:] for x in bulk.chrom_num]

    bulk["chrom"] = "chr" + bulk.chrom_num
    bulk.chrom_num = bulk.chrom_num.replace({"X": 23, "Y": 24})
    bulk.chrom_num = bulk.chrom_num.astype(int)
    bulk = bulk.sort_values(by=["chrom_num", "start"])
    bulk["CN"] = bulk["CN"] - 1
    bulk["CN"] = bulk["CN"].replace(0, 1)

    bulk["window"] = [f"{tup.chrom}:{tup.start}-{tup.end}" for tup in bulk.itertuples()]
    bulk = bulk.loc[bulk["chrom"] != "chrY"]
    bulk.set_index("window", inplace=True)
    return bulk


def double_figure(res, x, y1, y2):
    fig, ax = plt.subplots()
    # make a plot
    ax.plot(res[x], res[y1], color="red", marker="o")
    # set x-axis label
    ax.set_xlabel(x, fontsize=14)
    # set y-axis label
    ax.set_ylabel(y1, color="red", fontsize=14)

    ax2 = ax.twinx()
    # make a plot with different y-axis using second axis object
    ax2.plot(res[x], res[y2], color="blue", marker="o")
    ax2.set_ylabel(y2, color="blue", fontsize=14)
    plt.show()
