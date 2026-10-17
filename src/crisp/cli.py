"""Command-line interface for CRISP analyses and native calculations."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "native":
        from crisp.commands import main as native_main

        original = sys.argv
        try:
            sys.argv = [original[0], *original[2:]]
            native_main()
        finally:
            sys.argv = original
        return

    from crisp.analyses.runner import COMPONENTS, reproduce, run_component

    parser = argparse.ArgumentParser(
        prog="crisp",
        description="Cancer-gene Ranking Interpretation with Screen Profiles",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    core = commands.add_parser("reproduce", help="Derive the principal allocation summaries")
    core.add_argument("--data", type=Path, required=True, help="Extracted Supplementary_Data directory")
    core.add_argument("--output", type=Path, required=True, help="New output directory")
    core.add_argument("--union-readouts", action="store_true", help="Include continuous and label-sensitivity summaries")
    component = commands.add_parser("analysis", help="Derive one analysis or the manuscript figures")
    component.add_argument("component", choices=COMPONENTS)
    component.add_argument("--data", type=Path, required=True, help="Extracted Supplementary_Data directory")
    component.add_argument("--output", type=Path, required=True, help="New output directory")
    commands.add_parser("native", help="Native calculations; use 'crisp native --help' for commands")
    arguments = parser.parse_args()
    if arguments.command == "reproduce":
        reproduce(arguments.data, arguments.output, union_readouts=arguments.union_readouts)
    else:
        run_component(arguments.component, arguments.data, arguments.output)


if __name__ == "__main__":
    main()
