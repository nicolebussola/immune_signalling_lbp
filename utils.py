import numpy as np
import pandas as pd
from bokeh.layouts import column, layout, row
from bokeh.models import (BoxZoomTool, CDSView, ColorBar, ColumnDataSource,
                          CustomJS, CustomJSFilter, LassoSelectTool,
                          LinearColorMapper, RangeSlider, ResetTool, Slider,
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


def QC_metrics_UMAP_plot(adata):
    embedding = np.array(adata.obsm[f"X_umap"].astype(float))

    samples = np.array(list(adata.obs["log1p_total_counts"].index))
    field_counts = adata.obs["log1p_total_counts"].values.astype(float)
    field_mt = adata.obs["pct_counts_mt"].values.astype(float)
    field_ribo = adata.obs["pct_counts_ribo"].values.astype(float)
    field_scdblfinder_score = adata.obs["scDblFinder_score"].values.astype(float)
    field_scrublet_score = adata.obs["doublet_scores_scrublet"].values.astype(float)

    colormapper = LinearColorMapper(
        palette=Viridis256, low=min(field_counts), high=max(field_counts)
    )

    smp = np.expand_dims(samples, axis=1)
    data = np.hstack(
        (
            embedding,
            smp,
            np.expand_dims(field_counts, axis=1),
            np.expand_dims(field_mt, axis=1),
            np.expand_dims(field_ribo, axis=1),
            np.expand_dims(field_scdblfinder_score, axis=1),
            np.expand_dims(field_scrublet_score, axis=1),
        )
    )

    df = pd.DataFrame(
        data,
        columns=[
            "x",
            "y",
            "desc",
            "value_counts",
            "value_pct_mt",
            "value_pct_ribo",
            "scdblfinder_score",
            "scrublet_score",
        ],
    )

    df["x"] = pd.to_numeric(df["x"])
    df["y"] = pd.to_numeric(df["y"])

    p = figure(width=900, height=900, tooltips=TOOLTIPS, toolbar_location="left")

    p.title.align = "center"
    p.title.text_color = "black"
    p.title.text_font_size = "25px"

    p.add_tools(LassoSelectTool())
    p.add_tools(WheelZoomTool())
    p.add_tools(ZoomInTool())
    p.add_tools(ResetTool())
    p.add_tools(BoxZoomTool())

    p.hover.tooltips = """
        Sample name <strong>@desc</strong> <br>
        Index: $index <br>
        Coordinates: ($x, $y)<br>
        Counts <font face="Arial" size="2">@value_counts{0.2f} </font> <br> 
        Percent mito <font face="Arial" size="2">@value_pct_mt{0.2f} </font> <br> 
        Percent ribo <font face="Arial" size="2">@value_pct_ribo{0.2f} </font> <br> 
        scdbl score <font face="Arial" size="2">@scdblfinder_score{0.2f} </font> <br> 
        scrublet score <font face="Arial" size="2">@scrublet_score{0.2f} </font> <br> <br>
    
    """

    source = ColumnDataSource(data=df)

    slider_counts = RangeSlider(
        start=min(field_counts),
        end=max(field_counts),
        value=(min(field_counts), max(field_counts)),
        step=0.1,
        title="log1p total counts",
    )
    slider_mt = RangeSlider(
        start=min(field_mt),
        end=max(field_mt),
        value=(min(field_mt), max(field_mt)),
        step=1,
        title="% counts mito",
    )
    slider_ribo = RangeSlider(
        start=min(field_ribo),
        end=max(field_ribo),
        value=(min(field_ribo), max(field_ribo)),
        step=1,
        title="% counts ribo",
    )
    slider_scdbl = RangeSlider(
        start=min(field_scdblfinder_score),
        end=max(field_scdblfinder_score),
        value=(min(field_scdblfinder_score), max(field_scdblfinder_score)),
        step=0.1,
        title="scdblFinder score",
    )
    slider_scrublet = RangeSlider(
        start=min(field_scrublet_score),
        end=max(field_scrublet_score),
        value=(min(field_scrublet_score), max(field_scrublet_score)),
        step=0.1,
        title="scrublet score",
    )

    callback = CustomJS(
        args=dict(s=source),
        code="""
        s.change.emit();
    """,
    )

    slider_counts.js_on_change("value", callback)
    slider_mt.js_on_change("value", callback)
    slider_ribo.js_on_change("value", callback)
    slider_scdbl.js_on_change("value", callback)
    slider_scrublet.js_on_change("value", callback)

    filt_counts = CustomJSFilter(
        args=dict(slider=slider_counts),
        code="""
            var indices = [];
            var start = slider.value[0];
            var end = slider.value[1];
            for (var i=0; i < source.get_length(); i++){
                if (source.data['value_counts'][i] >= start && source.data['value_counts'][i] <= end){
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }
            return indices;
            """,
    )

    filt_mt = CustomJSFilter(
        args=dict(slider=slider_mt),
        code="""
            var indices = [];
            var start = slider.value[0];
            var end = slider.value[1];
    
            for (var i=0; i < source.get_length(); i++){
                if (source.data['value_pct_mt'][i] >= start && source.data['value_pct_mt'][i] <= end){
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }
            return indices;
            """,
    )

    filt_ribo = CustomJSFilter(
        args=dict(slider=slider_ribo),
        code="""
            var indices = [];
            var start = slider.value[0];
            var end = slider.value[1];
    
            for (var i=0; i < source.get_length(); i++){
                if (source.data['value_pct_ribo'][i] >= start && source.data['value_pct_ribo'][i] <= end){
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }
            return indices;
            """,
    )

    filt_scdbl = CustomJSFilter(
        args=dict(slider=slider_scdbl),
        code="""
            var indices = [];
            var start = slider.value[0];
            var end = slider.value[1];
    
            for (var i=0; i < source.get_length(); i++){
                if (source.data['scdblfinder_score'][i] >= start && source.data['scdblfinder_score'][i] <= end){
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }
            return indices;
            """,
    )

    filt_scrublet = CustomJSFilter(
        args=dict(slider=slider_scrublet),
        code="""
            var indices = [];
            var start = slider.value[0];
            var end = slider.value[1];
    
            for (var i=0; i < source.get_length(); i++){
                if (source.data['scrublet_score'][i] >= start && source.data['scrublet_score'][i] <= end){
                    indices.push(true);
                } else {
                    indices.push(false);
                }
            }
            return indices;
            """,
    )

    view = CDSView(
        source=source,
        filters=[filt_counts, filt_mt, filt_ribo, filt_scdbl, filt_scrublet],
    )
    # view_mt = CDSView(source=source, filters=[filt_mt])

    p.scatter(
        x="x",
        y="y",
        source=source,
        size=3,
        color={"field": "value_counts", "transform": colormapper},
        alpha=0.8,
        view=view,
    )

    cb = ColorBar(color_mapper=colormapper, location=(5, 6))

    # Adding the color bar to the right side
    p.add_layout(cb, "right")

    return(
        row(
            p,
            column(
                slider_counts, slider_mt, slider_ribo, slider_scdbl, slider_scrublet
            ),
        )
    )


def interactive_embedding_blood_brain(adata, LABEL, embedding_method="umap"):
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

    p = figure(width=900, height=900, tooltips=TOOLTIPS, toolbar_location="left")

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
                size=3,
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
            size=3,
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


def interactive_embedding(adata, LABEL, embedding_method="umap"):
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
    embedding = np.array(adata.obsm[f"X_{embedding_method}"].astype(float))

    p = figure(width=900, height=900, tooltips=TOOLTIPS, toolbar_location="left")

    p.title.align = "center"
    p.title.text_color = "black"
    p.title.text_font_size = "25px"

    p.add_tools(LassoSelectTool())
    p.add_tools(WheelZoomTool())
    p.add_tools(ZoomInTool())
    p.add_tools(ResetTool())
    p.add_tools(BoxZoomTool())

    # categorical label
    if f"{LABEL}_colors" in list(adata.uns.keys()):
        mycols = adata.uns[f"{LABEL}_colors"]
        myclasses = pd.unique(adata.obs[LABEL])

        for col, theclass in zip(mycols, myclasses):
            idx = np.where(np.array(list(adata.obs[LABEL])) == str(theclass))[
                0
            ].tolist()

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
            p.scatter(
                x="x",
                y="y",
                source=source,
                size=3,
                color=col,
                alpha=0.8,
                legend_label=theclass,
            )
            p.legend.location = "top_left"
            p.legend.click_policy = "hide"

        return p

    # continuous label
    else:
        p.hover.tooltips = """
            Sample name <strong>@desc</strong> <br>
            Index: $index <br>
            Coordinates: ($x, $y)<br>
            Target value <font face="Arial" size="2">@value{0.2f} </font> <br> <br>
        """
        field = adata.obs[LABEL].values.astype(float)
        colormapper = LinearColorMapper(
            palette=Viridis256, low=min(field), high=max(field)
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
            step=0.1,
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
        # slider.js_link('value', r.glyph, 'color')
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

        # Adding the color bar to the right side
        p.add_layout(cb, "right")
        lay = layout([slider, p])
        return lay
