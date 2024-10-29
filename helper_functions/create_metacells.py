import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances


def check_consistancy(new_df, cell_dict):
    """
    Check if all cells are assigned to metacells
    new_df: dataframe with metacells
    cell_dict: dictionary with cell to metacell assignment
    """
    consistent = True
    for mc in new_df.cl:
        cells = [x for x in cell_dict if cell_dict[x] == mc]
        if len(cells) != new_df[new_df.cl == mc].n_cells.values[0]:
            print("Incorrect assignment")
            consistent = False
            break
    return consistent


def find_match(clust_df, d_tresh, saturated=False):
    """
    Find cells to merge based on the distance between them
    clust_df: dataframe with metacells
    d_tresh: distance threshold
    saturated: if True, allow merging cells with sufficient coverage
    return: list of tuples with cells to merge and list of cells to remove
    """
    to_merge_df = clust_df.copy()
    pcs = [x for x in clust_df.columns if x.startswith("PC")]
    dist = pairwise_distances(to_merge_df[pcs], n_jobs=20)
    np.fill_diagonal(dist, 1e8)
    dist_tup = []
    remove = set()
    for idx1 in range(dist.shape[0]):
        idx2 = np.argmin(dist[idx1, :])
        d = dist[idx1, idx2]
        if (
            d < d_tresh
            and clust_df.iloc[idx1].pass_min == False
            and clust_df.iloc[idx2].pass_min == False
        ):
            # merge 2 cells with unsufficient coverage
            dist_tup.append((to_merge_df.index[idx1], to_merge_df.index[idx2], d))
        elif (
            d < d_tresh
            and clust_df.iloc[idx1].pass_min == False
            and clust_df.iloc[idx2].pass_max
        ):
            # merge 1 cell with unsufficient coverage with another with sufficient
            dist_tup.append((to_merge_df.index[idx1], to_merge_df.index[idx2], d))
        elif d >= d_tresh and clust_df.iloc[idx1].pass_min == False:
            # cell has small coverage but no other cells near to merge
            remove.add(to_merge_df.index[idx1])
        elif saturated and d < d_tresh and clust_df.iloc[idx1].pass_min == False:
            dist_tup.append((to_merge_df.index[idx1], to_merge_df.index[idx2], d))

    dist_tup.sort(key=lambda a: a[2])
    merged = set()
    filtered_tup = []
    for tup in dist_tup:
        if tup[0] in merged or tup[1] in merged:
            pass
        else:
            merged.add(tup[0])
            merged.add(tup[1])
            filtered_tup.append(tup)
    return (filtered_tup, remove)


def merge_cells(d_tresh, clust_df, key, adata, min_t=200, max_t=1100):
    """
    Merge cells with small coverage into metacells
    d_tresh: distance threshold
    clust_df: dataframe with metacells
    key: key in adata.obs with cell clusters
    adata: AnnData object
    min_t: minimal coverage to pass
    max_t: maximal coverage to pass
    return: dataframe with metacells
    """
    not_passed = clust_df.shape[0] - clust_df.pass_min.sum()
    print(f"metacells to process: {not_passed}")
    clust_df0 = clust_df.copy()
    results_dfs = [clust_df0]
    it = 0
    cell_dict = adata.obs[key].to_dict()
    stagnation = False
    pcs = [x for x in clust_df.columns if x.startswith("PC")]

    while not_passed != 0:
        print(f"Iteration {it}")
        # if stagnation happens allow to merge small cells with large ones
        filtered_tup, remove = find_match(clust_df0, d_tresh, stagnation)
        new_df = []
        max_id = clust_df0.cl.astype(int).max()
        merged = set()
        for tup in filtered_tup:
            tup1 = clust_df0.loc[tup[0]]
            tup2 = clust_df0.loc[tup[1]]
            merged.add(tup[0])
            merged.add(tup[1])

            new_id = str(max_id + 1)
            max_id = max_id + 1
            new_cov = tup1.coverage + tup2.coverage
            new_n_cells = tup1.n_cells + tup2.n_cells
            new_x = (tup1.umap_x + tup2.umap_x) / 2
            new_y = (tup1.umap_y + tup2.umap_y) / 2
            new_pcs = (
                np.array([getattr(tup1, x) for x in pcs])
                + np.array([getattr(tup2, x) for x in pcs]) * 0.5
            )
            new_t_min = new_cov > min_t
            new_t_max = new_cov <= max_t
            new_df.append(
                (
                    new_id,
                    new_cov,
                    new_n_cells,
                    new_x,
                    new_y,
                    *new_pcs,
                    new_t_min,
                    new_t_max,
                )
            )
            merged.add(tup1.cl)
            merged.add(tup2.cl)
            cells1 = [x for x in cell_dict if cell_dict[x] == str(tup1.cl)]
            cells2 = [x for x in cell_dict if cell_dict[x] == str(tup2.cl)]
            for cell in cells1 + cells2:
                cell_dict[cell] = new_id

        for tup in clust_df0.itertuples():
            if tup[0] not in merged and tup[0] not in remove:
                new_df.append(tup[1:])
        new_df = pd.DataFrame(new_df, columns=clust_df.columns)
        new_df.cl = new_df.cl.astype(str)
        results_dfs.append(new_df)
        not_passed_new = new_df.shape[0] - new_df.pass_min.sum()
        stagnation = not_passed_new == not_passed
        not_passed = not_passed_new
        print(f"metacells to process: {not_passed}")
        if stagnation:
            print("stagnation")
        mcs = new_df.cl.unique()
        for cell in cell_dict:
            if cell_dict[cell] not in mcs:
                cell_dict[cell] = "-1"

        it += 1
        clust_df0 = new_df
        assert check_consistancy(new_df, cell_dict)

    print(f"DONE with {new_df.shape[0]} metacells in total")
    return (results_dfs[-1], cell_dict)
