"""Portable full numerical reproduction without original-machine path dependencies."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    args = parser.parse_args()
    data = args.data.resolve()
    source_root = Path(__file__).resolve().parents[5]
    output = args.output.resolve()
    if output.exists() or output.is_relative_to(data) or output.is_relative_to(ROOT):
        raise ValueError("Use a new output directory; stored results are never overwritten")
    manifest = json.loads((data / "BUNDLE_MANIFEST.json").read_text(encoding="utf-8"))
    for record in manifest["files"]:
        path = data / record["path"]
        assert path.is_file() and sha(path) == record["sha256"], record["path"]
    commands = [
        ["analyse_union.py", "--archive", str(data / "inputs/archive"), "--crisp-source", str(source_root),
         "--reference", str(data / "inputs/reference"), "--output", str(output)],
        ["verify_union.py", "--reference", str(data / "inputs/reference"), "--results", str(output)]]
    for command in commands:
        command[0] = str(ROOT / command[0])
        subprocess.run([sys.executable, "-B", *command], check=True)
    exact = []
    for path in sorted((data / "results").iterdir()):
        if path.suffix == ".csv":
            assert sha(path) == sha(output / path.name), "Numerical mismatch: " + path.name
            exact.append(path.name)
    # Binomial-tail evaluation can differ in the last binary64 digits across runtimes.
    stored = json.loads((data / "results/family_details.json").read_text(encoding="utf-8"))
    rebuilt = json.loads((output / "family_details.json").read_text(encoding="utf-8"))
    coverage_checks = []
    for context in ("CRC", "ESCA"):
        lhs = stored[context]["q_monte_carlo_interval"].pop("binomial_coverage")
        rhs = rebuilt[context]["q_monte_carlo_interval"].pop("binomial_coverage")
        assert math.isclose(lhs, rhs, rel_tol=0, abs_tol=1e-14), context
        coverage_checks.append({"context": context, "stored": lhs, "reproduced": rhs,
                                "absolute_tolerance": 1e-14})
    assert stored == rebuilt, "Numerical mismatch: family_details.json"
    report = {"status": "PASS", "input_and_package_files_verified": len(manifest["files"]),
              "exact_numerical_files": exact, "comparison": "CSV: byte-identical SHA-256; family details: exact values except binomial-tail roundoff",
              "binomial_coverage_checks": coverage_checks, "native_fits": 0}
    (output / "portable_reproduction.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
