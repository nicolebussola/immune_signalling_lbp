import argparse
from pathlib import Path

from ..utils import readable_directory
from .run_preprocessing import run_preprocessing


def get_parser():
    """Returns the parser for all the args for run_ddqc function used via script/cli."""
    parser = argparse.ArgumentParser(prog="run_preprocessing")
    parser.add_argument(
        "--input-path",
        "-i",
        action=readable_directory,
        default=Path(),
        help="Path of the directory that contains cellranger outputs (h5) for each sample and ddqc metrics.",
    )
    parser.add_argument(
        "--output-path-plot",
        "-o",
        action=readable_directory,
        default=Path(),
        help="Output path: where to save UMAP embedding plot.",
    )
    parser.add_argument(
        "--tissue",
        "-t",
        choices=["blood", "brain"],
        type=str,
        default="blood",
        help="Type of tissue (brain or blood), default: blood.",
    )
    parser.add_argument(
        "--n-top-genes",
        "-g",
        type=int,
        default=2000,
        help="Number of genes for feature selection",
    )
    parser.add_argument(
        "--save-plots",
        "-s",
        type=bool,
        default=True,
        help="Save interactive QC plots",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    run_preprocessing(
        input_path=Path(args.input_path),
        output_path_plot=Path(args.output_path_plot),
        tissue=args.tissue,
        n_top_genes=args.n_top_genes,
        save_plots=args.save_plots,
    )
