import textwrap

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from gseapy import enrichment_map
from matplotlib.colors import Normalize
from bokeh.models import (
    CDSView,
    ColorBar,
    ColumnDataSource,
    CustomJS,
    CustomJSFilter,
    Label,
    LinearColorMapper,
    RangeSlider,
)
from bokeh.palettes import Viridis256
from bokeh.plotting import figure

TOOLTIPS = [
    ("index", "$index"),
    ("(x,y)", "($x, $y)"),
    ("desc", "@desc"),
]
TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,examine,help"


def interactive_embedding(
    adata,
    label,
    embedding_method="umap",
    width=900,
    height=900,
    labels_loc="outside",
    palette_cont=Viridis256,
    title_font_size="25px",
):
    """Plot interactive Bokeh embedding for scRNA data, coloured by label."""
    samples = np.array(list(adata.obs[label].index))
    embedding = np.array(adata.obsm[f"X_{embedding_method}"].astype(float))

    p = figure(
        title=label,
        width=width,
        height=height,
        tooltips=TOOLTIPS,
        tools=TOOLS,
        toolbar_location="left",
    )
    p.title.align = "center"
    p.title.text_color = "black"
    p.title.text_font_size = title_font_size

    if f"{label}_colors" in list(adata.uns.keys()):
        mycols = adata.uns[f"{label}_colors"]
        myclasses = adata.obs[label].unique().categories
        for col, theclass in zip(mycols, myclasses):
            idx = np.where(np.array(list(adata.obs[label])) == str(theclass))[0].tolist()
            smp = np.expand_dims(samples[idx], axis=1)
            data = np.hstack((embedding[idx,], smp))
            df = pd.DataFrame(data, columns=["x", "y", "desc"])
            df["x"] = pd.to_numeric(df["x"])
            df["y"] = pd.to_numeric(df["y"])
            df["class"] = theclass
            source = ColumnDataSource.from_df(data=df)
            p.hover.tooltips = """
            Sample name <strong>@desc</strong> <br>
            Index: $index <br>
            Coordinates: ($x, $y)<br>
            Class <font face="Arial" size="2">@class </font> <br> <br>
            """
            p.scatter(x="x", y="y", source=source, size=3, color=col, alpha=0.8, legend_label=theclass)
            p.legend.click_policy = "hide"
            legend = p.legend[0]
            if labels_loc == "on_data":
                centroid_x = np.mean(df["x"], axis=0)
                centroid_y = np.mean(df["y"], axis=0)
                labels = Label(
                    x=centroid_x,
                    y=centroid_y,
                    text=theclass,
                    border_line_color="black",
                    border_line_alpha=1.0,
                    background_fill_color="white",
                    background_fill_alpha=1.0,
                )
                labels.text_font_style = "bold"
                p.add_layout(labels)
        p.add_layout(legend, "right")
        return p

    else:
        p.hover.tooltips = """
            Sample name <strong>@desc</strong> <br>
            Index: $index <br>
            Coordinates: ($x, $y)<br>
            Target value <font face="Arial" size="2">@value{0.2f} </font> <br> <br>
        """
        field = adata.obs[label].values.astype(float)
        colormapper = LinearColorMapper(palette=palette_cont, low=min(field), high=max(field))
        smp = np.expand_dims(samples, axis=1)
        data = np.hstack((embedding, smp, np.expand_dims(field, axis=1)))
        df = pd.DataFrame(data, columns=["x", "y", "desc", "value"])
        df["x"] = pd.to_numeric(df["x"])
        df["y"] = pd.to_numeric(df["y"])
        source = ColumnDataSource(data=df)

        slider = RangeSlider(
            start=min(field), end=max(field), value=(min(field), max(field)), step=0.05, title="Value"
        )
        callback = CustomJS(args=dict(s=source), code="s.change.emit();")
        slider.js_on_change("value", callback)

        filt = CustomJSFilter(
            args=dict(slider=slider),
            code="""
                var indices = [];
                var start = slider.value[0];
                var end = slider.value[1];
                for (var i=0; i < source.get_length(); i++){
                    if (source.data['value'][i] >= start && source.data['value'][i] <= end){
                        indices.push(true);
                    } else {
                        indices.push(false);
                    }
                }
                return indices;
            """,
        )
        view = CDSView(source=source, filters=[filt])
        p.scatter(
            x="x", y="y", source=source, size=3,
            color={"field": "value", "transform": colormapper}, alpha=0.8, view=view,
        )
        cb = ColorBar(color_mapper=colormapper, location=(5, 6))
        p.add_layout(cb, "right")
        from bokeh.layouts import layout
        return layout([slider, p])


def _nudge_positions(pos, x_shift=0, y_shift=0.07):
    return {node: (x + x_shift, y + y_shift) for node, (x, y) in pos.items()}


def _wrap_text(text, width=30):
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False))


def _extract_ora_df(pdata, db_name, source_label):
    """
    Extract ORA estimates from pdata into a tidy DataFrame.

    Parameters:
    pdata (AnnData): Pseudo-bulk AnnData with ora_estimate_{db_name} in obsm.
    db_name (str): Database key.
    source_label (str): Label added as 'Source' column (e.g. 'Brain (Cohort 1)').

    Returns:
    pd.DataFrame with columns: pathway scores + 'cell_t' + 'Source'.
    """
    df = pdata.obsm[f"ora_estimate_{db_name}"].copy()
    df.index = df.index.str.replace("Astrocytes", "Astrocyte", regex=False)
    df["cell_t"] = [x.split("_")[1] if "_" in x else x for x in df.index]
    df = df[~df.index.str.contains("nan", na=False)]
    df["Source"] = source_label
    return df

def plot_pseudobulk_ora(sources, db_name, paths, out_folder):
    """
    Boxplot + swarmplot of pseudo-bulk ORA estimates for selected pathways.

    Parameters:
    sources (list of tuple): Each entry is (label, pdata) where label is a
        key in SOURCE_PALETTE (e.g. 'Brain (Cohort 1)').
    db_name (str): Database key ('wikipathways', 'reactome', or 'go').
    paths (list): Pathway names to plot.
    out_folder (Path): Directory to save figures.
    """
    dfs = [_extract_ora_df(pdata, db_name, label) for label, pdata in sources]
    combined_df = pd.concat(dfs, ignore_index=True)

    available = [p for p in paths if p in combined_df.columns]
    if not available:
        print(f"No matching pathways found in {db_name} ORA results — skipping plot.")
        return

    palette = {k: v for k, v in SOURCE_PALETTE.items() if k in combined_df["Source"].unique()}

    for pathway in available:
        fig, ax = plt.subplots(figsize=(max(8, combined_df["cell_t"].nunique() * 1.5), 4))

        sns.boxplot(
            data=combined_df, x="cell_t", y=pathway, hue="Source", palette=palette, ax=ax,
        )
        sns.swarmplot(
            data=combined_df, x="cell_t", y=pathway, hue="Source", palette=palette,
            size=3, dodge=True, edgecolor="#444445", linewidth=1, alpha=0.8,
            legend=False, ax=ax,
        )

        cell_types = combined_df["cell_t"].unique()
        for i in range(1, len(cell_types)):
            ax.axvline(i - 0.5, color="black", linestyle="--", linewidth=0.3)

        ax.set_title(pathway, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("ORA score")
        ax.tick_params(axis="x", rotation=0, labelsize=12)
        ax.legend(loc="upper right", fontsize=12, bbox_to_anchor=(1.2, 0.5))

        plt.tight_layout()
        safe_name = pathway[:60].replace("/", "_")
        fig.savefig(out_folder / f"pseudobulk_ora_{db_name}_{safe_name}.pdf", bbox_inches="tight")
        plt.close(fig)


def plot_prod_dict(prod_dict, out_file, figsize=(12, 9), cmap="BuPu"):
    """
    Plot LR product expression as a heatmap (cell types × LR pairs).

    Parameters:
    prod_dict (dict): {lr_pair: {cell_type: value}} from compute_prod_dict.
    out_file (str or Path): Output file path.
    figsize (tuple): Figure size. Default (12, 9).
    cmap (str): Colormap. Default "BuPu".
    """
    df = pd.DataFrame(prod_dict)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df, cmap=cmap, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=14)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=90, fontsize=18)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close(fig)
