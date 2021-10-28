from collections import defaultdict
from synister import SynisterDb
import argparse
import matplotlib.pyplot as plt
import pandas as pd


ALL_SPLITS = ["brain_region", "synapse", "skeleton", "hemi_lineage", "known"]


def plot_comparison(
    dataframe,
    splits,
    cmap=plt.get_cmap("tab20"),
    hatch="xx",
):
    PARTITION_ORDER = ["train", "test", "validation"]

    nplots = len(splits) + 1
    ncols = 3
    nrows = (nplots + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 4), sharey=True)
    for ax in axs.flat[nplots:]:
        ax.axis("off")

    datasets = dataframe["dataset"].unique()
    n_ds = len(datasets)
    nt_total = sorted(list(dataframe["nt"].dropna().unique()))
    n_nt_total = len(nt_total)
    nt_colors = {c: cmap(i) for i, c in enumerate(nt_total)}

    # Top-left plot is raw NT counts for each dataset.
    nt_ax = axs.flat[0]
    nt_counts = dataframe.groupby(["dataset", "nt"]).size()
    nt_counts.unstack(0).plot(kind="bar", ax=nt_ax)
    nt_ax.legend().set_visible(False)
    nt_ax.set(xlabel=None, xticks=[])
    h, _ = nt_ax.get_legend_handles_labels()

    for j, pa in enumerate(h):
        for i, rect in enumerate(pa.patches):
            rect.set_facecolor(cmap(i))
            rect.set_hatch(hatch * j)

    def order_sort(col, order):
        ordinal = pd.Categorical(col, categories=order, ordered=True)
        return pd.Series(ordinal)

    for ax, split in zip(axs.flat[1:], splits):
        split_col = f"split_{split}"
        if split_col not in dataframe.columns:
            print(f"Split `{split}` not in dataframe, skipping...")
            continue
        else:
            df = dataframe[dataframe[split_col].notna()]

        ds_n_nt = defaultdict(lambda: 0)
        for dataset, dfset in df.groupby(["dataset"]):
            nt_counts = dfset.groupby([split_col, "nt"]).size()
            nt_counts = nt_counts.unstack()
            nts = sorted(list(nt_counts.columns))
            ds_n_nt[dataset] = nt_counts.shape[1]
            nt_counts = nt_counts.sort_values(
                by=split_col, key=lambda col: order_sort(col, PARTITION_ORDER)
            )
            nt_counts.plot(
                kind="bar", stacked=True, ax=ax, color=[nt_colors[n] for n in nts]
            )

        ax.set_title(f"{split} split")
        ax.legend().set_visible(False)
        ax.set(xlabel=None)
        ax.tick_params(axis="x", labelrotation=45)

        # Group bars by dataset and style.
        h, _ = ax.get_legend_handles_labels()
        rect_width = 1 / float(n_ds + 1)
        for ds_i, n_nt in enumerate(ds_n_nt[n] for n in datasets):
            ds_h = h[:n_nt]
            h = h[n_nt:]
            for j, pa in enumerate(ds_h):
                for rect in pa.patches:
                    rect.set_x(
                        rect.get_x()
                        + ds_i * rect_width
                        - (rect.get_width() - rect_width) / 2.0
                    )
                    rect.set_hatch(hatch * ds_i)
                    rect.set_width(rect_width)

    fig.tight_layout()
    # Draw legends to side.
    (handles, labels) = axs.flat[1].get_legend_handles_labels()

    n = []
    for i in range(n_ds):
        n.append(axs.flat[-1].bar(0, 0, color="gray", hatch=hatch * i))

    ds_leg = plt.legend(n, dataframe["dataset"].unique(), loc=[1.05, 0.1])
    plt.legend(handles[:n_nt_total], labels[:n_nt_total], loc=[1.05, 0.4])
    fig.add_artist(ds_leg)

    return fig


def load_datasets(credentials, db_names):
    dbs = {db_name: SynisterDb(credentials, db_name) for db_name in db_names}

    records = {
        name: {
            "skeletons": db.get_skeletons(),
            "synapses": db.get_synapses(),  # No need to fetch splits, as this will fetch the union set.
        }
        for name, db in dbs.items()
    }

    return records


def datasets_to_dataframe(datasets, drop_missing_skids=True):
    dataframe = None
    for name, dataset in datasets.items():
        skel_df = pd.DataFrame.from_dict(dataset["skeletons"], orient="index")
        skel_df = skel_df.loc[
            skel_df.index.dropna()
        ]  # Must `dropna` to get rid of errant `None` skeleton IDs
        skel_df.index = skel_df.index.astype("int64")
        syn_df = pd.DataFrame.from_dict(dataset["synapses"], orient="index")
        syn_df["dataset"] = name

        if drop_missing_skids:
            syn_df = syn_df[syn_df["skeleton_id"].notna()]
        syn_df["nt"] = (
            skel_df.loc[syn_df["skeleton_id"], "nt_known"]
            .map(lambda x: x[0] if x is not None else None)
            .values
        )

        syn_df = pd.concat(
            [
                syn_df,
                pd.DataFrame(
                    syn_df["splits"].values.tolist(), index=syn_df.index
                ).add_prefix("split_"),
            ],
            axis=1,
        )

        if dataframe is None:
            dataframe = syn_df
        else:
            dataframe = dataframe.append(syn_df)

    return dataframe


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot a comparison of neurotransmitter counts between two datasets, broken out by split partitions"
    )
    parser.add_argument(
        "dataset_databases",
        type=str,
        nargs="+",
        help="MongoDB database names for each dataset to compare",
    )
    parser.add_argument("--credentials", "-c", type=str, help="MongoDB credential file")
    parser.add_argument(
        "--splits",
        "-s",
        type=lambda s: list(s.split(",")),
        help="Comma-separated list of split names to analyze",
    )

    args = parser.parse_args()

    records = load_datasets(args.credentials, args.dataset_databases)
    dataframe = datasets_to_dataframe(records)

    fig = plot_comparison(records, args.splits)

    filename = f"{'_'.join(args.dataset_databases.join)}_{'_'.join(args.splits)}.svg"
    fig.savefig(filename)
