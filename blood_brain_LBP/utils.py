import argparse
import os

import numpy as np
import pandas as pd
from bokeh.layouts import column, grid, layout, row
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
from bokeh.plotting import figure, output_file, show

from .labels import QC_LABELS_SAMPLE

TOOLTIPS = [
    ("index", "$index"),
    ("(x,y)", "($x, $y)"),
    ("desc", "@desc"),
]
TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,examine,help"


class readable_directory(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        _dir = values
        if not os.path.isdir(_dir):
            raise argparse.ArgumentTypeError(f"{_dir} is not a valid path.")
        if os.access(_dir, os.R_OK):
            setattr(namespace, self.dest, _dir)
        else:
            raise argparse.ArgumentTypeError(f"Directory: {_dir} is not a readable")


def detect_column_type(column, threshold=10):
    unique_values = column.unique()
    num_unique = len(unique_values)

    if num_unique <= threshold:
        return "Discrete"
    return "Numeric"


def QC_metrics_UMAP_plot(adata):
    embedding = np.array(adata.obsm[f"X_umap"].astype(float))

    samples = np.array(list(adata.obs.index))
    extracted_fields = {}
    metrics = QC_LABELS_SAMPLE + [c for c in adata.obs.columns if "doublet" in c]
    for label in metrics:
        adata_col = adata.obs[label]
        col_type = detect_column_type(adata_col)
        if col_type == "Numeric":
            extracted_fields[label] = adata_col.values.astype(float)
        else:
            extracted_fields[label] = adata_col.values.astype("str")

    data_arrays = (
        [embedding]
        + [np.expand_dims(samples, axis=1)]
        + [np.expand_dims(extracted_fields[field], axis=1) for field in metrics]
    )

    data = np.hstack(data_arrays)

    df = pd.DataFrame(
        data,
        columns=["x", "y", "sample"] + metrics,
    )

    df["x"] = pd.to_numeric(df["x"])
    df["y"] = pd.to_numeric(df["y"])
    source = ColumnDataSource(data=df)
    callback = CustomJS(
        args=dict(s=source),
        code="""
        s.change.emit();
    """,
    )
    p = figure(
        width=900, height=900, tools=TOOLS, tooltips=TOOLTIPS, toolbar_location="left"
    )

    p.title.align = "center"
    p.title.text_color = "black"
    p.title.text_font_size = "25px"

    p.hover.tooltips = """
        Sample name <strong>@desc</strong> <br>
        Index: $index <br>
        Coordinates: ($x, $y)<br>
        Counts <font face="Arial" size="2">@value_counts{0.2f} </font> <br>
        Percent mito <font face="Arial" size="2">@value_pct_mt{0.2f} </font> <br>
        Percent ribo <font face="Arial" size="2">@value_pct_ribo{0.2f} </font> <br>
        Scdbl score <font face="Arial" size="2">@scdblfinder_score{0.2f} </font> <br>
        Scrublet score <font face="Arial" size="2">@scrublet_score{0.2f} </font> <br> 
        Number of genes with positive counts <font face="Arial" size="2">@log1p_n_genes_by_counts{0.2f} </font> <br>
        Cumulative percentage of counts for 20 most expressed genes <font face="Arial" size="2">@pct_counts_in_top_20_genes{0.2f} </font> <br> <br>
  
    """
    sliders = []
    QC_LABELS_CONT = [
        label for label in metrics if detect_column_type(adata.obs[label]) == "Numeric"
    ]
    for col_name in QC_LABELS_CONT:
        min_val = min(df[col_name].astype(float))
        max_val = max(df[col_name].astype(float))
        slider = RangeSlider(
            start=min_val,
            end=max_val,
            value=(min_val, max_val),
            step=0.1,
            title=col_name,
        )
        slider.js_on_change("value", callback)
        sliders.append(slider)

    callback = CustomJS(
        args=dict(s=source),
        code="""
        s.change.emit();
    """,
    )

    filt_list = [
        CustomJSFilter(
            args=dict(slider=slider),
            code=f"""
                        var indices = [];
                        var start = slider.value[0];
                        var end = slider.value[1];
                        for (var i=0; i < source.get_length(); i++){{
                            if (source.data['{col_name}'][i] >= start && source.data['{col_name}'][i] <= end){{
                                indices.push(true);
                            }} else {{
                                indices.push(false);
                            }}
                        }}
                        return indices;
                    """,
        )
        for col_name, slider in zip(QC_LABELS_CONT, sliders)
    ]

    view = CDSView(
        source=source,
        filters=filt_list,
    )

    colormapper = LinearColorMapper(
        palette=Viridis256,
        low=min(df[QC_LABELS_CONT[0]].astype(float)),
        high=max(df[QC_LABELS_CONT[0]].astype(float)),
    )
    p.scatter(
        x="x",
        y="y",
        source=source,
        size=3,
        color={"field": QC_LABELS_CONT[0], "transform": colormapper},
        alpha=0.8,
        view=view,
    )

    cb = ColorBar(color_mapper=colormapper, location=(5, 6))
    p.add_layout(cb, "right")

    return row(
        p,
        column(*sliders),
    )


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
    """plot interactive plot for scRNA data

    Args:
        adata (AnnData): scanpy AnnData object. Must have the "tissue" (blood/brain) in adata.obs
        label (str): target to color the points
        embedding_method (string): name of embedding vectors in adata.obsm, e.g. UMAP


    Returns:
        bokeh.plotting._figure.figure: interactive embedding plot colored by label
    """
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

    # categorical label
    if f"{label}_colors" in list(adata.uns.keys()):
        mycols = adata.uns[f"{label}_colors"]
        myclasses = adata.obs[label].unique().categories
        for col, theclass in zip(mycols, myclasses):
            idx = np.where(np.array(list(adata.obs[label])) == str(theclass))[
                0
            ].tolist()

            smp = np.expand_dims(samples[idx], axis=1)

            data = np.hstack(
                (
                    embedding[idx,],
                    smp,
                )
            )
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
            p.scatter(
                x="x",
                y="y",
                source=source,
                size=3,
                color=col,
                alpha=0.8,
                legend_label=theclass,
            )
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

    # continuous label
    else:
        p.hover.tooltips = """
            Sample name <strong>@desc</strong> <br>
            Index: $index <br>
            Coordinates: ($x, $y)<br>
            Target value <font face="Arial" size="2">@value{0.2f} </font> <br> <br>
        """
        field = adata.obs[label].values.astype(float)
        colormapper = LinearColorMapper(
            palette=palette_cont, low=min(field), high=max(field)
        )
        smp = np.expand_dims(samples, axis=1)
        data = np.hstack((embedding, smp, np.expand_dims(field, axis=1)))
        df = pd.DataFrame(data, columns=["x", "y", "desc", "value"])
        df["x"] = pd.to_numeric(df["x"])
        df["y"] = pd.to_numeric(df["y"])
        source = ColumnDataSource(data=df)

        slider = RangeSlider(
            start=min(field),
            end=max(field),
            value=(min(field), max(field)),
            step=0.05,
            title="Value",
        )
        callback = CustomJS(
            args=dict(s=source),
            code="""
            s.change.emit();
        """,
        )
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
            x="x",
            y="y",
            source=source,
            size=3,
            color={"field": "value", "transform": colormapper},
            alpha=0.8,
            view=view,
        )

        cb = ColorBar(color_mapper=colormapper, location=(5, 6))
        p.add_layout(cb, "right")
        lay = layout([slider, p])
        return lay


def gridlayout(
    labels, adata, width=800, height=700, ncols=3, fname=None, labels_loc="outside"
):
    ps = []
    for label in labels:
        p = interactive_embedding(
            adata, label, width=width, height=height, labels_loc=labels_loc
        )
        ps.append(p)

    if fname is not None:
        output_file(fname)

    show(grid(ps, ncols=ncols))
