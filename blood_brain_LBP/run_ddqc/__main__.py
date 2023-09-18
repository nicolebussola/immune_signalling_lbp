import argparse
from pathlib import Path

from ..utils import readable_directory
# from .run_ddqc import run_ddqc


def get_parser():
    """Returns the parser for all the args for run_ddqc function use via script/cli."""
    parser = argparse.ArgumentParser(prog="run_ddqc")
    parser.add_argument(
        "--input-path",
        "-i",
        action=readable_directory,
        default=Path(),
        help="Input path of the directory that contains ...",  # TODO: expand help
    )
    parser.add_argument(
        "--output-path",
        "-o",
        action=readable_directory,
        default=Path(),
        help="Output path: where to store the results.",
    )
    parser.add_argument(
        "--tissue",
        "-t",
        choices=["blood",],
        type=str,
        default="blood",
        help="Type of tissue, default: blood.",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    # run_ddqc(
    #     input_path=Path(args.input_path),
    #     output_path=Path(args.output_path),
    #     tissue=args.tissue,
    # )
