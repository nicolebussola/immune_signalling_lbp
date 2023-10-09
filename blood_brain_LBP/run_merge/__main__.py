import argparse
from pathlib import Path

from ..utils import readable_directory
from .merge_data import run_merge


def get_parser():
    """Returns the parser for all the args for run_ddqc function used via script/cli."""
    parser = argparse.ArgumentParser(prog="run_merge")
    parser.add_argument(
        "--input-path",
        "-i",
        action=readable_directory,
        default=Path(),
        help="Path of the directory that contains h5ad for each sample after QC.",
    )
    parser.add_argument(
        "--batch",
        "-b",
        choices=["1", "2"],
        type=str,
        help="Batch corresponding to input data",
    )
    parser.add_argument(
        "--output-path",
        "-o",
        action=readable_directory,
        default=Path(),
        help="Output path: where to store the resulting merged data.",
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
        "--method-hvg",
        "-m",
        choices=["HighlyDeviant", "cell_ranger"],
        type=str,
        default="cell_ranger",
        help="Method for feature selection, default: cell_ranger.",
    )
    parser.add_argument(
        "--n-top-genes",
        "-n",
        type=int,
        default=4000,
        help="Number of genes for feature selection. Default 4000",
    )
    parser.add_argument(
        "--output-path-plot",
        "-p",
        action=readable_directory,
        default=Path(),
        help="Where to save the plots.",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    run_merge(
        input_path=Path(args.input_path),
        batch=args.batch,
        output_path=Path(args.output_path),
        tissue=args.tissue,
        method_hvg=args.method_hvg,
        n_top_genes=args.n_top_genes,
        output_path_plot=Path(args.output_path_plot),
    )
