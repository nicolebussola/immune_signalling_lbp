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
