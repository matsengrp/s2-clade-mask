#!/usr/bin/env python3
"""Pipeline rule entrypoint: run one ``matUtils extract`` pass on the MAT.

Wraps ``s2clademask.run_sample_paths`` / ``run_clade_paths`` so the documented
``matUtils`` ``-d /`` output-dir quirk (see the comment block in the library) is
never dropped by inlining the command into the Snakefile shell. Emits either the
root->tip sample-paths TSV (``--kind sample``; large, tens of GB) or the
per-clade founder clade-paths TSV (``--kind clade``).

    python scripts/extract_paths.py --config config.yaml --kind sample --out results/sample_paths.tsv
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml",
                   help="Path to the YAML config (default: config.yaml).")
    p.add_argument("--kind", choices=["sample", "clade"], required=True,
                   help="Which matUtils extract pass to run.")
    p.add_argument("--out", required=True, help="Output TSV path.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    log_file = cfg["log"]
    tree = cfg["tree"]  # required by load_config

    out_dir = os.path.dirname(os.path.abspath(cli.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if cli.kind == "sample":
        lib.run_sample_paths(tree, cli.out, cfg["matutils"], log_file)
    else:
        lib.run_clade_paths(tree, cli.out, cfg["matutils"], log_file)
    lib.log(f"Wrote {cli.kind}-paths to {cli.out}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
