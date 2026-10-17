"""Command-line access to the CRISP supplementary-data calculations."""
import argparse
from pathlib import Path

from crisp.analyses.runner import reproduce


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--union-readouts", action="store_true")
    args = parser.parse_args()
    reproduce(args.data, args.output, args.union_readouts)


if __name__ == "__main__":
    main()
