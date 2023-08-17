import os
from pathlib import Path

import ddqc
import pandas as pd
import pegasus as pg
import pegasusio as io

BATCH = 1
tissue = "blood"


if BATCH == 1:
    batch_path = Path("../LBP_brain_blood_pairs/data/narsad_cellRanger_outs/")
else:
    batch_path = Path()


for PT in [
    "PT-182",
    "PT-185",
    "PT-201",
    "PT-203",
    "PT-205",
    "PT-206",
    "PT-208",
    "PT-212",
    "PT-214",
]:
    for side in ["R", "L"]:
        try:
            if tissue == "blood":
                h5_name = f"{PT}-{side}-B_CellBender_filtered.h5"
            else:
                h5_name = f"{PT}-{side}_CellBender_filtered.h5"

            print("\n======================")
            adata = io.read_input(
                batch_path / tissue / f"{PT}-{tissue}-{side}" / h5_name, genome="hg38"
            )
            print(adata.shape)
            data_qc = ddqc.ddqc_metrics(
                adata,
                return_df_qc=True,
                display_plots=False,
                ribo_prefix="^RP[SL]|^RPLP|^RPSA",
            )
            pg.filter_data(adata)

            print(adata.shape)
            #h5_out_name = os.path.splitext(h5_name)[0] + "_ddqc.h5"
            # pg.write_output(
            #     adata, batch_path / tissue / f"{PT}-{tissue}-{side}" / h5_out_name
            # )
            data_qc.to_csv(
                batch_path
                / tissue
                / f"{PT}-{tissue}-{side}"
                / f"{PT}-{side}-B_CellBender_filtered_ddqc.csv"
            )
            print("======================\n")
        except Exception as e:
            print(e)
