#!/usr/bin/env python3
"""Pipeline rule entrypoint: compare the two anchoring modes and emit the final YAML.

Reads the Mode 1 (representative) and Mode 2 (founder) node-map artifacts produced
by ``resolve_nodes.py`` and does two things:

1. **Comparison (report-only).** Builds the full Nextstrain-clade mask map under
   each anchoring and diffs them (``diff_masks``), plus a per-PANGO-clade node
   comparison, writing it all to ``--out-comparison``. This never fails the run --
   per design invariant #1 the authoritative output stays representative-primary;
   divergence is surfaced for human review, not treated as an error.

2. **Final output.** Resolves the authoritative node map as **Mode 1 with a
   per-clade Mode-2 fallback only where the config approves it**
   (``apply_fallbacks``; ``representative_fallbacks: {via: founder|sample}``). An
   unresolved clade with no approved fallback is a hard error. Builds the masks and
   writes ``--out-yaml`` (``usher_masked_sites.yaml``).

    python scripts/compare_and_build.py -c config.yaml \\
        --nodes-representative results/nodes_representative.yaml \\
        --nodes-founder results/nodes_founder.yaml \\
        --sample-paths results/sample_paths.tsv --clade-paths results/clade_paths.tsv \\
        --out-yaml results/usher_masked_sites.yaml \\
        --out-comparison results/anchoring_comparison.txt
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def write_comparison_report(path, spec, rep_nodes, rep_via, rep_unresolved,
                            founder_nodes, founder_unresolved, rep_masks,
                            founder_masks):
    """Write the human-readable rep-vs-founder comparison report (report-only)."""
    lines = []
    lines.append("=" * 78)
    lines.append("ANCHORING COMPARISON  (representative / Mode 1  vs  founder / Mode 2)")
    lines.append("Report only -- the authoritative output is representative-primary "
                 "(invariant #1).")
    lines.append("=" * 78)

    # --- Node-level: where do the two anchorings differ per PANGO clade? -----
    lines.append("\n## Per-PANGO-clade anchor node")
    same, diff = [], []
    for clade in sorted(spec):
        r = rep_nodes.get(clade)
        f = founder_nodes.get(clade)
        if r is not None and f is not None and r == f:
            same.append(clade)
        else:
            diff.append(clade)
            note = ""
            if r is None:
                note += "  (rep unresolved)"
            if f is None:
                note += "  (no founder)"
            lines.append(f"  {clade}: rep={r}  founder={f}{note}")
    lines.append(f"  -> {len(same)} clade(s) resolve to the SAME node either way; "
                 f"{len(diff)} differ (or one side is unavailable).")
    if rep_unresolved:
        lines.append(f"  Mode 1 could not resolve: {sorted(rep_unresolved)}")
    if founder_unresolved:
        lines.append(f"  Mode 2 has no founder for: {sorted(founder_unresolved)}")

    # --- Mask-level: do the emitted Nextstrain-clade masks differ? -----------
    lines.append("\n## Emitted Nextstrain-clade masks")
    diffs = lib.diff_masks(rep_masks, founder_masks)
    if not diffs:
        lines.append("  IDENTICAL: every Nextstrain clade gets the same in-range "
                     "sites and mutations under both anchorings.")
    else:
        lines.append(f"  DIFFERENT: {len(diffs)} Nextstrain clade(s) differ:")
        for clade, e in diffs.items():
            lines.append(f"\n  {clade}:")
            if e["sites_only_rep"]:
                lines.append(f"    sites only under REP:     {sorted(e['sites_only_rep'])}")
            if e["sites_only_founder"]:
                lines.append(f"    sites only under FOUNDER: {sorted(e['sites_only_founder'])}")
            if e["muts_only_rep"]:
                lines.append(f"    mutations only under REP:     {sorted(e['muts_only_rep'])}")
            if e["muts_only_founder"]:
                lines.append(f"    mutations only under FOUNDER: {sorted(e['muts_only_founder'])}")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml")
    p.add_argument("--nodes-representative", help="Mode 1 node-map YAML.")
    p.add_argument("--nodes-founder", help="Mode 2 node-map YAML.")
    p.add_argument("--sample-paths", help="matUtils --sample-paths TSV.")
    p.add_argument("--clade-paths", help="matUtils --clade-paths TSV.")
    p.add_argument("--out-yaml", help="Final usher_masked_sites.yaml.")
    p.add_argument("--out-comparison", help="Comparison report text file.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    log_file = cfg["log"]
    lib.validate_resolution_overrides(cfg["representative_name_overrides"],
                                      cfg["representative_fallbacks"])

    nodes_rep = cli.nodes_representative or cfg["nodes_representative"]
    nodes_founder = cli.nodes_founder or cfg["nodes_founder"]
    sample_paths = cli.sample_paths or cfg["sample_paths"]
    clade_paths = cli.clade_paths or cfg["clade_paths"]
    out_yaml = cli.out_yaml or cfg["output"]
    out_comparison = cli.out_comparison or cfg["comparison_report"]
    if not sample_paths:
        raise ValueError("need --sample-paths (or config 'sample_paths')")
    if not clade_paths:
        raise ValueError("need --clade-paths (or config 'clade_paths')")

    for path in (out_yaml, out_comparison):
        d = os.path.dirname(os.path.abspath(path))
        if d:
            os.makedirs(d, exist_ok=True)

    spec, version = lib.parse_branch_spec(cfg["branch_spec"])
    lib.log(f"Parsed {len(spec)} PANGO clades, spec version {version}", log_file)
    clade_founders = lib.parse_clade_paths(clade_paths)
    lib.log(f"Parsed {len(clade_founders):,} clade annotations from {clade_paths}",
            log_file)
    ns_full_values = lib.load_nextstrain_clade_values(
        cfg["metadata"], cfg["clade_column"], log_file)
    start, end = cfg["site_include_start"], cfg["site_include_end"]

    rep_nodes, rep_via, rep_unresolved = lib.read_nodes(nodes_rep)
    founder_nodes, _, founder_unresolved = lib.read_nodes(nodes_founder)
    lib.log(f"Loaded Mode 1 ({len(rep_nodes)} nodes, {len(rep_unresolved)} "
            f"unresolved) from {nodes_rep}", log_file)
    lib.log(f"Loaded Mode 2 ({len(founder_nodes)} nodes, {len(founder_unresolved)} "
            f"unresolved) from {nodes_founder}", log_file)

    # --- Comparison (report-only) --------------------------------------------
    lib.log("\n=== Building masks under representative anchoring (Mode 1) ===", log_file)
    rep_masks = lib.build_masks_for_nodes(
        spec, rep_nodes, ns_full_values, clade_founders, sample_paths, start, end, log_file)
    lib.log("\n=== Building masks under founder anchoring (Mode 2) ===", log_file)
    founder_masks = lib.build_masks_for_nodes(
        spec, founder_nodes, ns_full_values, clade_founders, sample_paths, start, end, log_file)
    write_comparison_report(out_comparison, spec, rep_nodes, rep_via, rep_unresolved,
                            founder_nodes, founder_unresolved, rep_masks, founder_masks)
    lib.log(f"Wrote anchoring comparison to {out_comparison}", log_file)

    # --- Final output: Mode 1 + approved per-clade fallback ------------------
    lib.log("\n=== Resolving final node map (representative primary; approved "
            "fallbacks only) ===", log_file)
    final_nodes = dict(rep_nodes)
    final_via = dict(rep_via)
    final_nodes, final_via = lib.apply_fallbacks(
        spec, final_nodes, final_via, clade_founders, sample_paths,
        cfg["representative_fallbacks"], log_file)

    lib.log("\n=== Building final masks ===", log_file)
    clade_masks = lib.build_masks_for_nodes(
        spec, final_nodes, ns_full_values, clade_founders, sample_paths, start, end, log_file)

    # Report contributing PANGO clades per Nextstrain clade.
    lib.log("Mask sources per Nextstrain clade (contributing PANGO clades):", log_file)
    for ns_clade in sorted(clade_masks):
        sources = [s for s in clade_masks[ns_clade]["sources"] if not s.startswith("old:")]
        if sources:
            lib.log(f"  {ns_clade} <- {sources}", log_file)

    # union_old force-include must happen *before* writing (mutates clade_masks);
    # the strict recovery gate itself lives in the separate check_recovery rule.
    if cfg["union_old"] and cfg["old_yaml"]:
        old_sites = lib.load_old_sites(cfg["old_yaml"])
        lib.check_recovers_old(old_sites, clade_masks, union_old=True,
                               include_range=(start, end), log_file=log_file)

    dropped_exclusions = {c: e["exclusions"] for c, e in spec.items() if e["exclusions"]}
    provenance = [
        f"source spec:   {cfg['branch_spec']}",
        f"tree:          {os.path.basename(cfg['tree']) if cfg['tree'] else '(sample-paths supplied)'}",
        f"metadata:      {os.path.basename(cfg['metadata'])}",
        f"clade column:  {cfg['clade_column']}",
        f"include range: [{start}, {end}] (masks outside this range are omitted)",
    ]
    lib.write_output(out_yaml, clade_masks, spec, provenance, dropped_exclusions)
    lib.log(f"Wrote {out_yaml}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
