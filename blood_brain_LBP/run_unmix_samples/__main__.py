import argparse
from pathlib import Path

from ..utils import readable_directory
from .unmix_samples import run_unmix_samples


def get_parser():
    """Returns the parser for all the args for run_ddqc function used via script/cli."""
    parser = argparse.ArgumentParser(prog="run_ddqc")
    parser.add_argument(
        "--project-path",
        "-i",
        action=readable_directory,
        default=Path(),
        help="Path of the directory that contains the .h5ad merged, QC-ed data of blood - batch 2",
    )
    parser.add_argument(
        "--blood-adata-name-mixed",
        "-f",
        type=str,
        default="",
        help="Name of the merged blood data file",
    )

    parser.add_argument(
        "--blood-adata-name",
        "-fo",
        type=str,
        default="",
        help="Name of the data file where to saved unmerged data",
    )
    parser.add_argument(
        "--output-path-plot",
        "-o",
        action=readable_directory,
        default=Path(),
        help="Output path: where to save UMAP plot.",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    run_unmix_samples(
        project_path=Path(args.project_path),
        blood_adata_name_mixed=args.blood_adata_name_mixed,
        output_path_plot=args.output_path_plot,
    )
