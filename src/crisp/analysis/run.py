import csv
import gzip
import json
from fractions import Fraction
from pathlib import Path

from . import cross_method, dawn_factorial, domain_bridge, factorial_decisions, prodigy_factorial

def object_pairs(pairs):
    output = {}
    for (key, value) in (pairs):
        if (key in output):
            raise ValueError("Duplicate JSON key: " + key)
        output[key] = value
    return output

def reject_constant(value):
    raise ValueError("Nonfinite JSON token: " + value)

def read_input(base, item):
    specification = {"path": item, "format": "json"} if (isinstance(item, str)) else item
    if ("paths" in specification):
        return [read_input(base, value) for (value) in (specification["paths"])]
    if ("mapping" in specification):
        return {key: read_input(base, value) for (key, value) in (specification["mapping"].items())}
    path = (base / specification["path"]).resolve()
    opener = gzip.open if (path.suffix == ".gz") else open
    with opener(path, "rt", encoding = "utf-8-sig") as stream:
        kind = specification.get("format", "json")
        if (kind == "text"):
            return stream.read()
        if (kind in ("csv", "tsv")):
            return list(csv.DictReader(stream, delimiter = "\t" if (kind == "tsv") else ","))
        if (kind == "jsonl"):
            return [
                json.loads(line, object_pairs_hook = object_pairs, parse_constant = reject_constant)
                for (line) in (stream)
                if (line.strip())
            ]
        if (kind != "json"):
            raise ValueError("Unsupported declared input format")
        result = json.load(
            stream, object_pairs_hook = object_pairs, parse_constant = reject_constant
        )
        for (key) in (specification.get("select", [])):
            result = result[key]
        return result

def serializable(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if (isinstance(value, dict)):
        if (any(not isinstance(key, str) for (key) in (value))):
            raise ValueError("Output requires explicit string keys or row records")
        return {key: serializable(item) for ((key, item)) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [serializable(item) for (item) in (value)]
    if (isinstance(value, set)):
        return [serializable(item) for (item) in (sorted(value))]
    return value

def run_operation(config_path, registry):
    config_path = Path(config_path).resolve()
    config = read_input(config_path.parent, config_path.name)
    if (config["operation"] not in registry):
        raise ValueError("Unknown scientific operation")
    inputs = config.get("inputs", {})
    options = config.get("options", {})
    if (set(inputs) & set(options)):
        raise ValueError("An argument cannot be both a file input and an option")
    arguments = {
        name: read_input(config_path.parent, item) for ((name, item)) in (inputs.items())
    }
    arguments.update(options)
    output = (config_path.parent / config["output"]).resolve()
    if (output.exists()):
        raise FileExistsError("Use a distinct output path; existing results are retained")
    result = registry[config["operation"]](**arguments)
    output.parent.mkdir(parents = True, exist_ok = True)
    with output.open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(serializable(result), stream, ensure_ascii = False, allow_nan = False, indent = 2)
        stream.write("\n")
    return {"operation": config["operation"], "output": str(output)}

OPERATIONS = {
    "prodigy_freeze": prodigy_factorial.freeze,
    "prodigy_evaluate": prodigy_factorial.evaluate,
    "factorial_decision_freeze": factorial_decisions.freeze,
    "factorial_decision_evaluate": factorial_decisions.evaluate,
    "dawn_freeze": dawn_factorial.build_freeze,
    "dawn_evaluate": dawn_factorial.evaluate,
    "cross_method": cross_method.analysis,
    "bridge_freeze": domain_bridge.prepare_bridge,
    "bridge_evaluate": domain_bridge.finish_from_prepared,
}

def run(config_path):
    return run_operation(config_path, OPERATIONS)
