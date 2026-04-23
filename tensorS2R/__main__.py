import argparse
from pathlib import Path

import scanpy as sc
from .adata_processing_utils import prepare_scran_de_data_micro_mono
from .factorization_run import run_tensorcell2cell
from .enrichment_run import (
    create_enrichment_directory,
    enrich_mono_micro_blood_brain,
    enrich_mono_micro_blood_micro,
    enrich_mono_micro_brain,
)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--project_path",
        help="Path were data are stored",
        default="/sc/arion/projects/psychgen/lbp/data/ScProcesses_brainBlood_nicole/",
    )
    parser.add_argument(
        "-b", "--cohort", help="cohort identifier", choices=["cohort_1", "cohort_2"]
    )
    parser.add_argument(
        "-s", "--sample_key", help="Patient column in adata.obs", default="sample"
    )
    parser.add_argument(
        "-g",
        "--groupby",
        help="Celltype column in adata.obs for liana",
        default="cell_type_coarse",
    )
    parser.add_argument(
        "-f",
        "--factorization_type",
        help="Strategy for tensor factorization. This name is used to create out folder",
        choices=["brain_coarse", "micro_blood_coarse", "brain_blood_coarse"],
    )
    parser.add_argument(
        "-m", "--merging_mode", help="Merging mode for c2c", default="inner"
    )
    parser.add_argument(
        "-r", "--factorization_rank", help="c2c factorization rank", type=int, default=None
    )
    parser.add_argument("-d", "--device", help="run on GPU or CPU", default="cuda")
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    project_path = Path(args.project_path)

    try:
        adata_br = sc.read(
            project_path
            / args.cohort
            / "c2c_liana_outputs"
            / "data_brain_scran_de_micro_mono.h5ad"
        )
        adata_bl = sc.read(
            project_path
            / args.cohort
            / "c2c_liana_outputs"
            / "data_blood_scran_de_micro_mono.h5ad"
        )
    except FileNotFoundError:
        adata_bl, adata_br = prepare_scran_de_data_micro_mono(
            args.cohort, project_path, "cell_type.v2", "cell_type.v2"
        )
    loadings_folder, adata_tf = run_tensorcell2cell(
        project_path,
        args.cohort,
        args.sample_key,
        args.groupby,
        args.factorization_type,
        args.merging_mode,
        args.factorization_rank,
        args.device,
    )

    enrich_path = create_enrichment_directory(loadings_folder)
    if args.factorization_type == "brain_coarse":
        enrich_mono_micro_brain(
            args.cohort,
            loadings_folder,
            adata_br,
            genedbs=["Reactome_2022", "GO_Biological_Process_2023"],
            design=args.factorization_type,
            out_folder=enrich_path,
        )

    elif args.factorization_type == "micro_blood_coarse":
        ## UPLOAD ALREADY DEG COMPUTED
        adata_tf = sc.read(
            project_path / args.cohort / "c2c_liana_outputs" / "blood_microglia_tf.h5ad"
        )
        enrich_mono_micro_blood_micro(
            args.cohort,
            loadings_folder,
            adata_br,
            adata_bl,
            adata_tf,
            design=args.factorization_type,
            out_folder=enrich_path,
        )

    elif args.factorization_type == "brain_blood_coarse":
        ## UPLOAD ALREADY DEG COMPUTED
        adata_tf = sc.read(
            project_path / args.cohort / "c2c_liana_outputs" / "blood_brain_tf.h5ad"
        )
        enrich_mono_micro_blood_brain(
            args.cohort,
            loadings_folder,
            adata_br,
            adata_bl,
            adata_tf,
            design=args.factorization_type,
            out_folder=enrich_path,
        )
