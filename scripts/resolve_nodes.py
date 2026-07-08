#!/usr/bin/env python3
"""Pipeline rule entrypoint: resolve PANGO clades to tree nodes in one anchoring mode.

Two modes, written as separate node-map artifacts (``write_nodes`` schema) so the
downstream ``compare_and_build`` rule can diff them and build the final output:

* ``--mode representative`` (Mode 1, string match): anchor each PANGO clade at its
  representative sample's node via an exact pre-``|`` name match (applying the
  approved ``representative_name_overrides``, but **no** founder fallback). Clades
  whose representative is absent from the tree are recorded as ``unresolved``.
* ``--mode founder`` (Mode 2): anchor each PANGO clade at its annotated founder
  (MRCA) node from ``--clade-paths``. Clades with no annotated founder are
  recorded as ``unresolved``.

    python scripts/resolve_nodes.py -c config.yaml --mode representative \\
        --sample-paths results/sample_paths.tsv --out results/nodes_representative.yaml
    python scripts/resolve_nodes.py -c config.yaml --mode founder \\
        --clade-paths results/clade_paths.tsv --out results/nodes_founder.yaml
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
    p.add_argument("--mode", choices=["representative", "founder"], required=True,
                   help="Anchoring mode to resolve.")
    p.add_argument("--sample-paths", help="matUtils --sample-paths TSV (mode representative).")
    p.add_argument("--clade-paths", help="matUtils --clade-paths TSV (mode founder).")
    p.add_argument("--out", required=True, help="Output node-map YAML path.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    log_file = cfg["log"]
    lib.validate_resolution_overrides(cfg["representative_name_overrides"],
                                      cfg["representative_fallbacks"])
    spec, version = lib.parse_branch_spec(cfg["branch_spec"])
    lib.log(f"Parsed {len(spec)} PANGO clades, spec version {version}", log_file)

    out_dir = os.path.dirname(os.path.abspath(cli.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if cli.mode == "representative":
        sample_paths = cli.sample_paths or cfg["sample_paths"]
        if not sample_paths:
            raise ValueError("mode representative needs --sample-paths (or config 'sample_paths')")
        lib.log("Mode 1 (representative / string match): resolving representative "
                "sample nodes...", log_file)
        clade_to_node, via, unresolved = lib.resolve_representative_only(
            spec, sample_paths, cfg["representative_name_overrides"], log_file)
        meta = {"mode": "representative", "spec_version": version,
                "source": os.path.abspath(sample_paths)}
    else:
        clade_paths = cli.clade_paths or cfg["clade_paths"]
        if not clade_paths:
            raise ValueError("mode founder needs --clade-paths (or config 'clade_paths')")
        clade_founders = lib.parse_clade_paths(clade_paths)
        lib.log(f"Parsed {len(clade_founders):,} clade annotations from {clade_paths}",
                log_file)
        lib.log("Mode 2 (founder / MRCA): resolving annotated founder nodes...", log_file)
        clade_to_node, no_founder = lib.resolve_founder_nodes(spec, clade_founders)
        via = {c: "founder" for c in clade_to_node}
        unresolved = no_founder
        meta = {"mode": "founder", "spec_version": version,
                "source": os.path.abspath(clade_paths)}
        if no_founder:
            lib.log(f"  {len(no_founder)} clade(s) have no annotated founder node: "
                    f"{no_founder}", log_file)

    lib.write_nodes(cli.out, clade_to_node, via=via, unresolved=unresolved, meta=meta)
    lib.log(f"Resolved {len(clade_to_node)}/{len(spec)} clades ({len(unresolved)} "
            f"unresolved) -> {cli.out}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
