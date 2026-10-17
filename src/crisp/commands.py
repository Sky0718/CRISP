"""Native calculations and exact allocation summaries for CRISP."""
import argparse
from importlib import import_module
from pathlib import Path

COMMANDS = {
    "prepare": ("crisp.core.prepare", "prepare_common_domain"),
    "prepare-molecular": ("crisp.core.molecular_run", "prepare_molecular"),
    "compute": ("crisp.core.run", "run_job"),
    "export-native": ("crisp.core.export", "export_results"),
    "analyse": ("crisp.analysis.run", "run"),
    "precision": ("crisp.precision.workflow", "run"),
    "verify-precision": ("crisp.precision.verification", "run"),
    "mechanism": ("crisp.mechanism.workflow", "run"),
    "oesophageal": ("crisp.oesophageal.workflow", "run"),
    "statistics": ("crisp.statistics.workflow", "run"),
    "diagnostics": ("crisp.diagnostics.workflow", "run"),
}

def main():
    parser = argparse.ArgumentParser(prog="crisp.py native")
    parser.add_argument("command", choices=tuple(COMMANDS))
    parser.add_argument("configuration", type=Path)
    arguments = parser.parse_args()
    module, function = COMMANDS[arguments.command]
    return getattr(import_module(module), function)(arguments.configuration.resolve())
