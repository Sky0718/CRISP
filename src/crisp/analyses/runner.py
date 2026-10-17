"""Portable entry points for CRISP analyses with external supplementary data.

Input data are read-only. Each command requires a new output directory.
Scientific validation outputs are written beside newly derived results.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
SOURCE_ROOT = PACKAGE_ROOT.parent
CHECKOUT_ROOT = (
    SOURCE_ROOT.parent
    if SOURCE_ROOT.name == "src" and (SOURCE_ROOT.parent / "pyproject.toml").is_file()
    else PACKAGE_ROOT
)
COMPONENTS = (
    "union", "union-ledger", "continuous", "categories", "composition",
    "candidates", "structure", "opportunity", "alleles", "graph-endpoints",
    "graph-numerical", "figures",
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepare(data: str | Path, output: str | Path) -> tuple[Path, Path]:
    source, out = Path(data).resolve(), Path(output).resolve()
    if not source.is_dir() or not (source / "Union_Sensitivity_Analysis").is_dir():
        raise ValueError("Data must be the extracted Supplementary Data directory.")
    if out.exists() or out.is_relative_to(source) or out.is_relative_to(CHECKOUT_ROOT):
        raise ValueError("Choose a new output directory outside the input and code directories.")
    out.mkdir(parents=True)
    return source, out


def _python(script: str, *arguments: str | Path) -> None:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [sys.executable, "-B", str(HERE / script), *map(str, arguments)],
        check=True, env=environment,
    )


def _check_input_map(base: Path) -> int:
    with (base / "INPUT_OUTPUT_MAP.csv").open(encoding="utf-8-sig", newline="") as stream:
        records = list(csv.DictReader(stream))
    for record in records:
        relative = Path(record["package_path"])
        path = (base / relative).resolve()
        if not path.is_relative_to(base.resolve()):
            raise ValueError("Input manifest path escapes its data directory.")
        if _sha(path) != record["sha256"]:
            raise ValueError("Input checksum mismatch: " + relative.as_posix())
    return len(records)


def _compare_csvs(expected: Path, actual: Path) -> list[str]:
    names = []
    for path in sorted(expected.glob("*.csv")):
        if path.read_bytes() != (actual / path.name).read_bytes():
            raise AssertionError("Numerical table differs: " + path.name)
        names.append(path.name)
    if not names:
        raise ValueError("Expected CSV reference tables are missing.")
    return names


def _continuous(source: Path, out: Path) -> None:
    _python("Union_Validation/reproduce.py", "--data", source / "Union_Validation",
            "--work", out,
            "--archive", source / "Reader_Data/Data_10/Union_Analysis/inputs/archive")


def reproduce(data: str | Path, output: str | Path, union_readouts: bool = False) -> Path:
    """Reconstruct cancer-listed binary summaries and check figure-source hashes.

    ``union_readouts=True`` also derives the continuous and label-sensitivity
    outputs. The separately supplied reference tables remain unchanged.
    """
    source, out = _prepare(data, output)
    extension = source / "Union_Sensitivity_Analysis"
    input_count = _check_input_map(extension)
    _python("Union_Sensitivity_Analysis/analyse_fraction.py", "--union",
            extension / "inputs/union", "--output", out / "binary_union")
    tables = _compare_csvs(extension / "expected/results", out / "binary_union")
    mapping = json.loads((source / "FIGURE_SOURCE_MAP.json").read_text(encoding="utf-8"))
    for item in mapping["records"]:
        path = (source / item["file"]).resolve()
        if not path.is_relative_to(source) or _sha(path) != item["packaged_sha256"]:
            raise AssertionError("Figure-source checksum mismatch: " + item["file"])
    if union_readouts:
        _continuous(source, out / "continuous_union")
    report = {
        "input_files": input_count, "binary_tables": tables,
        "figure_source_records": len(mapping["records"]),
        "continuous_union": union_readouts,
    }
    (out / "checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return out


def run_component(component: str, data: str | Path, output: str | Path) -> Path:
    """Derive one named analysis without writing into the supplied data tree."""
    if component not in COMPONENTS:
        raise ValueError("Unknown component. Choose: " + ", ".join(COMPONENTS))
    source, out = _prepare(data, output)
    archive = source / "Reader_Data/Data_10/Union_Analysis/inputs/archive"
    if component == "union":
        _python("Union_Sensitivity_Analysis/reproduce.py", "--data",
                source / "Union_Sensitivity_Analysis", "--output", out / "results")
    elif component == "union-ledger":
        folder = source / "Reader_Data/Data_10/Union_Analysis"
        _python("Reader_Data/Data_10/Union_Analysis/reproduce.py", "--data", folder,
                "--output", out / "results")
    elif component == "continuous":
        _continuous(source, out / "results")
    elif component == "categories":
        _python("Colorectal_Functional_Analysis/Code/derive_categories.py",
                "--archive", archive,
                "--crisp-source", HERE.parent.parent,
                "--reference-tables", source / "functional_content",
                "--input-references", source / "Colorectal_Functional_Analysis/Code/INPUT_REFERENCES.json",
                "--output", out / "results")
    elif component == "composition":
        _python("Allocation_and_Graph_Controls/Composition/Code/derive_composition.py",
                "--source", source, "--output", out / "results")
        _python("Allocation_and_Graph_Controls/Composition/Code/independent_check.py",
                "--source", source, "--output", out / "results")
    elif component == "candidates":
        folder = source / "Allocation_and_Graph_Controls/Candidates"
        _python("Allocation_and_Graph_Controls/Candidates/Code/reproduce_candidate_context.py",
                "--source-data", folder, "--output", out / "results")
        _python("Allocation_and_Graph_Controls/Candidates/Code/audit_candidate_context.py",
                "--source-data", folder, "--report", out / "candidate_checks.json")
    elif component == "structure":
        folder = source / "Allocation_and_Graph_Controls/Graph_Controls"
        result = out / "results"
        result.mkdir()
        _python("Variant_and_Network_Context/Structure/derive_structure.py",
                "--source", folder, "--output", result)
        _python("Variant_and_Network_Context/Structure/verify_structure.py",
                "--source", folder, "--results", result, "--report", out / "structure_checks.json")
    elif component == "opportunity":
        _python("Opportunity_and_Oesophageal_Analysis/Structure/check_opportunity_denominators.py",
                "--source", source / "Variant_and_Network_Context/Structure/query_gene_records.json",
                "--output", out / "opportunity_checks.json")
    elif component == "alleles":
        # Cached public source records and exact identifiers accompany the data.
        # Work on a private copy because the source-specific script writes beside them.
        folder = source / "Variant_and_Network_Context/Alleles"
        work = out / "alleles"
        shutil.copytree(folder, work, ignore=shutil.ignore_patterns("*.py", "__pycache__"))
        for name in ("annotate_coverage.py", "verify_coverage.py", "PROTOCOL.txt"):
            shutil.copy2(HERE / "Variant_and_Network_Context/Alleles" / name, work)
        _python("Variant_and_Network_Context/Alleles/annotate_coverage.py",
                "summarize", "--output-dir", work)
        _python("Variant_and_Network_Context/Alleles/verify_coverage.py", "--work-dir", work)
    elif component == "graph-endpoints":
        _python("Allocation_and_Graph_Controls/Graph_Controls/rederive_endpoints.py",
                "--directory", source / "Allocation_and_Graph_Controls/Graph_Controls",
                "--output", out / "endpoint_checks.json")
    elif component == "graph-numerical":
        _python("Allocation_and_Graph_Controls/Graph_Controls/verify_numerical_portable.py",
                "--directory", source / "Allocation_and_Graph_Controls/Graph_Controls",
                "--output", out / "numerical_checks.json")
    elif component == "figures":
        from .figure_panels import run as render_source_panels

        _python("Main_Figure_Reproduction/build_figures.py", "--inputs",
                source / "Main_Figure_Reproduction/inputs", "--output", out / "figures")
        work = out / "composition"
        work.mkdir()
        _python("Variant_and_Network_Context/Figures/build_composition_figure.py",
                "--source", source / "Allocation_and_Graph_Controls/Composition",
                "--output", out / "figures/Fig_Composition.pdf", "--work", work)
        render_source_panels(source, out / "figures")
    return out
