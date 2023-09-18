import argparse
from pathlib import Path

from ..utils import readable_directory
from .run_ddqc import run_ddqc


def get_parser():
    """Returns the parser for all the args for run_ddqc function used via script/cli."""
    parser = argparse.ArgumentParser(prog="run_ddqc")
    parser.add_argument(
        "--input-path",
        "-i",
        action=readable_directory,
        default=Path(),
        help="Path of the directory that contains cellranger outputs (h5) for each sample.",
    )
    parser.add_argument(
        "--output-path",
        "-o",
        action=readable_directory,
        default=Path(),
        help="Output path: where to store the resulting csv with ddqc metrics.",
    )
    parser.add_argument(
        "--tissue",
        "-t",
        choices=[
            "blood",
        ],
        type=str,
        default="blood",
        help="Type of tissue (brain or blood), default: blood.",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    run_ddqc(
        input_path=Path(args.input_path),
        output_path=Path(args.output_path),
        tissue=args.tissue,
    )
