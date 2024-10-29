"""Cell calling based on barcode rank plot."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm


def prepare_barcodes(cellranger_output: str) -> pd.DataFrame:
    """Prepare dataframe with fragments per barcode.

    Parameters
    ----------
    cellranger_output
        Path to '*fragments.tsv.gz' file in Cellranger's output.

    Returns
    -------
    Dataframe with fragments per barcode.
    """

    filename = cellranger_output
    print("test")

    cols_to_use = [3, 4]
    accumulated_counts = pd.Series(dtype="int")

    # Process in chunks
    chunksize = 10**6
    for chunk in tqdm(
        pd.read_csv(
            filename,
            chunksize=chunksize,
            comment="#",
            sep="\t",
            header=None,
            usecols=cols_to_use,
        )
    ):
        chunk_counts = chunk.groupby(3)[4].sum()
        accumulated_counts = accumulated_counts.add(chunk_counts, fill_value=0)

    accumulated_counts = pd.DataFrame(accumulated_counts, columns=["counts"])
    accumulated_counts.index.name = "barcodes"

    return accumulated_counts


def _barcode_rank_plot(sort_counts, i_max, figure_path):
    """Plot barcode rank plot."""

    hist, bins = np.histogram(sort_counts, bins=np.arange(1, np.max(sort_counts) + 2))
    cumulative_hist = np.cumsum(hist[::-1])[::-1]

    x1 = cumulative_hist[cumulative_hist <= i_max]
    y1 = bins[:-1][cumulative_hist <= i_max]

    x2 = cumulative_hist[cumulative_hist > i_max]
    y2 = bins[:-1][cumulative_hist > i_max]

    plt.figure(figsize=(10, 6))

    plt.scatter(x1, y1, color="blue", label=f"Cells")
    plt.scatter(x2, y2, color="grey", label=f"Background")

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Barcodes")
    plt.ylabel("Counts")
    plt.title("Barcode rank plot")
    plt.legend()
    plt.grid(True)
    plt.savefig(figure_path, dpi=300)


def compute_highest_gradient(
    count_matrix: pd.DataFrame,
    hi_limit: float = 99.99,
    lo_limit: int = 2000,
    figure_path: str = "barcode_rank_plot.png",
) -> pd.DataFrame:
    """Compute the highest gradient on the barcode rank plot.

    Parameters
    ----------
    counts
        Dataframe with barcodes and their fragment count.
    hi_limit
        The higher cutoff for barcodes took into consideration for gradient.
        By default, 99.99th percentile is used here.
    lo_limit
        The lower cutoff for barcodes took into consideration for gradient.
        By default, 2000 fragments is used here as cells under 2000 fragments
        will be useless for downstream analysis.
    figure_path
        Path to the directory where you want to store your barcode rank plot figure.

    Returns
    -------
    List with barcodes above the highest gradient.
    """

    sort_counts = np.array(sorted(count_matrix["counts"])[::-1])

    hi = np.percentile(sort_counts, hi_limit)
    lo = lo_limit

    filtered_counts = np.log10(sort_counts[(sort_counts < hi) & (sort_counts > lo)])

    i_max = np.argmax(np.abs(np.gradient(filtered_counts)))
    barcode_cutoff = len(sort_counts[sort_counts >= 10 ** (filtered_counts[i_max])])

    _barcode_rank_plot(
        sort_counts=sort_counts, i_max=barcode_cutoff, figure_path=figure_path
    )

    return count_matrix.sort_values(by="counts", ascending=False)["counts"][
        :barcode_cutoff
    ]
