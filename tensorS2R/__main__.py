import argparse
from pathlib import Path

import scanpy as sc

from .adata_processing_utils import prepare_scran_de_data_micro_mono
from .step1_factorize import run_tensorcell2cell
from .step2_enrich import run_enrichment


def get_parser():
    parser = argparse.ArgumentParser(
        description="tensorS2R pipeline — LIANA + Tensor-cell2cell, LR enrichment, and ORA"
    )
    parser.add_argument("-p", "--project_path", required=True, help="Root data directory")
    parser.add_argument("-b", "--cohort", required=True, choices=["cohort_1", "cohort_2"])
    parser.add_argument("-s", "--sample_key", default="sample", help="Patient column in adata.obs")
    parser.add_argument(
        "-g", "--groupby", default="cell_type_coarse", help="Cell-type column for LIANA predictions"
    )
    parser.add_argument(
        "-f", "--factorization_type",
        choices=["brain_coarse", "micro_blood_coarse", "brain_blood_coarse"],
        help="Analytical design (required for steps 1 and 2)",
    )
    parser.add_argument("-m", "--merging_mode", default="inner", help="Tensor merging mode")
    parser.add_argument(
        "-r", "--factorization_rank", type=int, default=None,
        help="Tensor rank; required when running step 2 without step 1",
    )
    parser.add_argument("-d", "--device", default="cuda", choices=["cuda", "cpu"])
    parser.add_argument(
        "-t", "--threshold_targets", type=float, default=0.1,
        help="Top fraction of genes for ORA (step 3; default: 0.1)",
    )
    parser.add_argument(
        "--step", default="all", choices=["1", "2", "3", "all"],
        help="Pipeline step to run: 1=factorize, 2=enrich, 3=ORA, all=run all (default: all)",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    project_path = Path(args.project_path)

    run_step1 = args.step in ("1", "all")
    run_step2 = args.step in ("2", "all")
    run_step3 = args.step in ("3", "all")

    if (run_step1 or run_step2) and args.factorization_type is None:
        parser.error("--factorization_type (-f) is required for steps 1 and 2.")
    if run_step2 and not run_step1 and args.factorization_rank is None:
        parser.error(
            "--factorization_rank (-r) is required when running step 2 without step 1."
        )

    # Load DE adatas needed for step 2
    if run_step2:
        cohort_liana_path = project_path / args.cohort / "c2c_liana_outputs"
        try:
            adata_br = sc.read(cohort_liana_path / "data_brain_scran_de_micro_mono.h5ad")
            adata_bl = sc.read(cohort_liana_path / "data_blood_scran_de_micro_mono.h5ad")
        except FileNotFoundError:
            adata_bl, adata_br = prepare_scran_de_data_micro_mono(
                args.cohort, project_path, "cell_type.v2", "cell_type.v2"
            )

        adata_tf = None
        if args.factorization_type == "micro_blood_coarse":
            adata_tf = sc.read(cohort_liana_path / "blood_microglia_tf.h5ad")
        elif args.factorization_type == "brain_blood_coarse":
            adata_tf = sc.read(cohort_liana_path / "blood_brain_tf.h5ad")

    # Step 1: Tensor factorization
    if run_step1:
        loadings_folder, _ = run_tensorcell2cell(
            project_path,
            args.cohort,
            args.sample_key,
            args.groupby,
            args.factorization_type,
            args.merging_mode,
            args.factorization_rank,
            args.device,
        )
    elif run_step2:
        # Reconstruct path from explicit rank when step 1 was already run separately
        loadings_folder = (
            project_path / args.cohort / "c2c_liana_outputs"
            / args.factorization_type
            / f"rank_{args.factorization_rank}"
            / args.merging_mode
            / args.groupby
        )

    # Step 2: LR filtering and enrichment
    if run_step2:
        run_enrichment(
            cohort=args.cohort,
            factorization_type=args.factorization_type,
            loadings_folder=loadings_folder,
            adata_br=adata_br,
            adata_bl=adata_bl,
            adata_tf=adata_tf,
        )

    # Step 3: ORA (lazy import to avoid loading MSigDB when not needed)
    if run_step3:
        from .step3_downstream import main as run_ora
        _lf = loadings_folder if (run_step1 or run_step2) else None
        run_ora(project_path, args.cohort, args.threshold_targets, _lf, args.factorization_type)
