import numpy as np
import pandas as pd
from bokeh.models import (BoxZoomTool, ColorBar, ColumnDataSource,
                          LassoSelectTool, LinearColorMapper, ResetTool,
                          WheelZoomTool, ZoomInTool)
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from bokeh.transform import factor_mark

TOOLTIPS = [
    ("index", "$index"),
    ("(x,y)", "($x, $y)"),
    ("desc", "@desc"),
]

MARKERS = ["hex", "x"]
TISSUES = ["brain", "blood"]

def interactive_embedding_blood_brain(adata, LABEL, embedding_method='umap'):
    """plot interactive plot for scRNA data

    Args:
        adata (AnnData): scanpy AnnData object. Must have the "tissue" (blood/brain) in adata.obs
        embedding_method (string): name of embedding vectors in adata.obsm, e.g. UMAP
        LABEL (str): target to color the points

    Returns:
        bokeh.plotting._figure.figure: interactive embedding plot colored by label
    """
    print(LABEL)
    samples = np.array(list(adata.obs[LABEL].index))
    tissues = np.array(list(adata.obs["tissue"].values))
    embedding = np.array(adata.obsm[f"X_{embedding_method}"].astype(float))


    p = figure(width=1200, height=1200, tooltips=TOOLTIPS, toolbar_location="left")

    p.title.align = "center"
    p.title.text_color = "black"
    p.title.text_font_size = "25px"

    # categorical label
    if f"{LABEL}_colors" in list(adata.uns.keys()):
        mycols = adata.uns[f"{LABEL}_colors"]
        myclasses = pd.unique(adata.obs[LABEL])

        for col, theclass in zip(mycols, myclasses):
            idx = np.where(np.array(list(adata.obs[LABEL])) == str(theclass))[
                0
            ].tolist()

            smp = np.expand_dims(samples[idx], axis=1)
            clas = np.expand_dims(tissues[idx], axis=1)

            data = np.hstack((embedding[idx,], smp, clas))
            df = pd.DataFrame(data, columns=["x", "y", "desc", "tissue"])
            df["x"] = pd.to_numeric(df["x"])
            df["y"] = pd.to_numeric(df["y"])
            source = ColumnDataSource.from_df(data=df)
            p.scatter(
                x="x",
                y="y",
                marker=factor_mark("tissue", MARKERS, TISSUES),
                source=source,
                size=5,
                color=col,
                alpha=0.8,
                legend_label=theclass,
            )
            p.legend.location = "top_left"
            p.legend.click_policy = "hide"

    # continuous label
    else:
        field = adata.obs[LABEL].values.astype(float)
        colormapper = LinearColorMapper(
            palette=Viridis256, low=min(field), high=max(field)
        )
        smp = np.expand_dims(samples, axis=1)
        clas = np.expand_dims(tissues, axis=1)
        data = np.hstack((embedding, smp, np.expand_dims(field, axis=1), clas))
        df = pd.DataFrame(data, columns=["x", "y", "desc", "target", "tissue"])
        df["x"] = pd.to_numeric(df["x"])
        df["y"] = pd.to_numeric(df["y"])
        source = ColumnDataSource.from_df(data=df)

        p.scatter(
            x="x",
            y="y",
            marker=factor_mark("tissue", MARKERS, TISSUES),

            source=source,
            size=5,
            color={"field": "target", "transform": colormapper},
            alpha=0.8,
        )
        cb = ColorBar(color_mapper=colormapper, location=(5, 6))

        # Adding the color bar to the right side
        p.add_layout(cb, "right")

    p.add_tools(LassoSelectTool())
    p.add_tools(WheelZoomTool())
    p.add_tools(ZoomInTool())
    p.add_tools(ResetTool())
    p.add_tools(BoxZoomTool())
    return p
