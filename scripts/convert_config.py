#!/usr/bin/env python3
"""Migrate the flat configuration preceding 15a2383 to the nested format."""

import argparse
from copy import deepcopy
from pathlib import Path
import sys

import yaml


KEY_PATHS = {
    "dir_run": ("run", "dir"),
    **{key: ("workflow", key) for key in ("mapper", "callers", "annotators")},
    **{key: ("inputs", key) for key in (
        "fasta", "index_minimap2", "dir_data", "suffix_fastq",
    )},
    **{key: ("filtering", key) for key in (
        "min_reads", "min_length_reads", "min_quality_mapping", "min_coverage",
        "min_size", "max_size", "min_dhffc", "max_dhbfc", "distance_sv",
        "n_callers", "consider_type", "consider_strand", "estimate_distance",
    )},
    **{key: ("annotation", key) for key in (
        "terms_relative", "species", "genome",
    )},
    **{f"{key}_{tool}": ("annotation", tool, key)
       for tool in ("snpeff", "vep", "annotsv")
       for key in ("version", "cache")},
    "max_size_vep": ("annotation", "vep", "max_size"),
    "bed_tandem_repeats": ("caller_settings", "tandem_repeats"),
    "bed_nvtr": ("caller_settings", "nvtr"),
    "config_nanosv": ("caller_settings", "nanosv", "config"),
    "bed_nanosv": ("caller_settings", "nanosv", "bed"),
    "model_clair3": ("caller_settings", "clair3", "model"),
    "model_svision": ("caller_settings", "svision", "model"),
    "resource_downloads": ("resources", "downloads"),
}


def convert_config(config):
    """Move legacy keys, retaining custom keys and existing nested settings.

    Missing keys stay missing. Mixed old/new input is rejected when both forms
    specify the same setting, so no value is silently discarded.
    """
    if not isinstance(config, dict):
        raise ValueError("configuration must be a YAML mapping")
    result = deepcopy({key: value for key, value in config.items()
                       if key not in KEY_PATHS})
    for old_key, path in KEY_PATHS.items():
        if old_key not in config:
            continue
        target = result
        for section in path[:-1]:
            target = target.setdefault(section, {})
            if not isinstance(target, dict):
                raise ValueError(f"cannot migrate {old_key}: {section} is not a mapping")
        if path[-1] in target:
            raise ValueError(f"both {old_key} and {'.'.join(path)} are present")
        target[path[-1]] = deepcopy(config[old_key])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="old YAML configuration")
    parser.add_argument("output", type=Path, nargs="?",
                        help="new file (must not exist); defaults to stdout")
    args = parser.parse_args()
    try:
        with args.input.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        converted = convert_config(config)
        text = yaml.safe_dump(converted, sort_keys=False, allow_unicode=True)
        if args.output is None:
            sys.stdout.write(text)
        else:
            with args.output.open("x", encoding="utf-8") as handle:
                handle.write(text)
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
