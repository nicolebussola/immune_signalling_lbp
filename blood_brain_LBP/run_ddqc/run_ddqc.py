import ddqc
import pegasusio as io
import logging

logging.basicConfig(format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S")


def run_ddqc(input_path, output_path, tissue):
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
                prefix = "B_" if tissue == "blood" else ""
                h5_name = f"{PT}-{side}-{prefix}CellBender_filtered.h5"
                adata = io.read_input(
                    input_path / tissue / f"{PT}-{tissue}-{side}" / h5_name,
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
                    / f"{PT}-{side}-{prefix}CellBender_filtered_ddqc.csv"
                )
            except Exception as e:
                logging.error(e)
