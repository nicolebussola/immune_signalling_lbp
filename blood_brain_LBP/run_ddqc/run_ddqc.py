import datetime
import logging
import os

import ddqc
import pandas as pd
import pegasusio as io

logging.getLogger().setLevel(logging.INFO)
logging.basicConfig(format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S")

timestamp = datetime.datetime.now().strftime("%m%d")


def run_ddqc(input_path, output_path, tissue):
    data_path = input_path / tissue
    patients = pd.unique(
        [
            "PT-" + p.split("-")[1]
            for p in os.listdir(data_path)
            if os.path.isdir(data_path / p)
        ]
    )
    logging.info(f"Patients found ({len(patients)}): {patients}")
    for PT in patients:
        for side in ["R", "L"]:
            try:
                prefix = "B_" if tissue == "blood" else ""
                h5_name = f"{PT}-{side}-{prefix}CellBender_filtered.h5"
                adata = io.read_input(
                    data_path / f"{PT}-{tissue}-{side}" / h5_name,
                    genome="hg38",
                )
                logging.info(f"adata shape before metrics: {adata.shape}")
                data_qc = ddqc.ddqc_metrics(
                    adata,
                    return_df_qc=True,
                    display_plots=False,
                    basic_n_genes=200,
                    ribo_prefix="^RP[SL]|^RPLP|^RPSA",
                )
                logging.info(f"adata shape after metrics: {adata.shape}")
                data_qc.to_csv(
                    output_path
                    / tissue
                    / f"{PT}-{tissue}-{side}"
                    / f"{PT}-{side}-{prefix}CellBender_filtered_ddqc_{timestamp}.csv"
                )
            except Exception as e:
                logging.error(e)
