import logging
import os

import ddqc
import pandas as pd
import pegasusio as io
from rich.logging import RichHandler

logging.getLogger().setLevel(logging.INFO)
FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


def run_ddqc(batch_path, tissue):
    data_path = batch_path / tissue
    patients = pd.unique(
        [
            "PT-" + p.split("-")[1]
            for p in os.listdir(data_path)
            if os.path.isdir(data_path / p)
        ]
    )
    log.info(f"Patients found ({len(patients)}): {patients}")
    for PT in sorted(patients):
        log.info(f"Computing metrics for patient: {PT}")
        for side in ["R", "L"]:
            patient_path = data_path / f"{PT}-{tissue}-{side}"
            log.info(f"Side {side}")
            try:
                h5_name = f"{PT}-{side}-CellBender_filtered.h5"
                adata = io.read_input(
                    patient_path / h5_name,
                    genome="hg38",
                )
                log.info(f"adata shape before metrics: {adata.shape}")
                log.info(f"Run ddqc")
                data_qc = ddqc.ddqc_metrics(
                    adata,
                    return_df_qc=True,
                    display_plots=False,
                    basic_n_genes=200,
                    ribo_prefix="^RP[SL]|^RPLP|^RPSA",
                )
                log.info(f"adata shape after metrics: {adata.shape}")
                data_qc.to_csv(
                    patient_path / f"{PT}-{side}_CellBender_filtered_ddqc.csv"
                )
            except Exception as e:
                log.exception(e)
