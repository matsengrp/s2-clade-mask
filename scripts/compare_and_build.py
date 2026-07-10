#!/usr/bin/env python3
"""Pipeline rule entrypoint: emit the founder-anchored masks + the anchor table.

Each PANGO clade is anchored at its annotated **founder (MRCA) node** (Mode 2). The
representative sample (Mode 1) is used only to *cross-check* that anchor. This rule:

1. **Builds two masking YAMLs** from the founder anchoring:
   * ``--out-yaml`` -- the mat-pcp **Nextstrain-clade** file (over-mask-coarsened),
   * ``--out-pango`` -- a companion **PANGO-clade** file (each clade accumulates its
     ancestral PANGO clades' masks; no Nextstrain coarsening).
   Both are built from a single shared sample-paths scan.

2. **Tabulates the founder-vs-representative distance** for every clade (patristic
   node/mutation distance + the exact set of mutations separating the two anchors)
   into ``--out-comparison`` (human) and ``--out-distances`` (machine-readable, read
   by the ``check_anchors`` rule which enforces ``max_anchor_distance``). This rule
   itself never fails on distance -- it always produces the table.

    python scripts/compare_and_build.py -c config.yaml \\
        --nodes-representative results/nodes_representative.yaml \\
        --nodes-founder results/nodes_founder.yaml \\
        --sample-paths results/sample_paths.tsv --clade-paths results/clade_paths.tsv \\
        --out-yaml results/usher_masked_sites.yaml \\
        --out-pango results/usher_masked_sites_pango.yaml \\
        --out-comparison results/anchoring_comparison.txt \\
        --out-distances results/anchor_distances.yaml
"""

import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def write_distance_report(path, distances, max_dist, absent_ok):
    """Human-readable table of founder-vs-representative anchor distances."""
    lines = []
    lines.append("=" * 78)
    lines.append("ANCHOR CROSS-CHECK  (founder anchor  vs  representative anchor)")
    lines.append("Output is founder-anchored; the representative is used only to verify")
    lines.append(f"the founder is within max_anchor_distance = {max_dist} mutations "
                 "(patristic).")
    lines.append("=" * 78)
    lines.append("")
    lines.append(f"{'clade':11} {'d_nodes':>7} {'d_muts':>6} {'status':7} relation / notes")
    worst = 0
    for clade in sorted(distances):
        d = distances[clade]
        if d.get("rep_absent"):
            ok = "OK" if clade in absent_ok else "FAIL"
            note = ("representative absent from tree (approved)" if clade in absent_ok
                    else "representative ABSENT and not in representative_absent_ok")
            lines.append(f"{clade:11} {'-':>7} {'-':>6} {ok:7} {note}")
            continue
        if "error" in d:
            lines.append(f"{clade:11} {'?':>7} {'?':>6} {'FAIL':7} {d['error']}")
            continue
        dm = d["d_muts"]
        worst = max(worst, dm)
        status = "OK" if dm <= max_dist else "FAIL"
        lines.append(f"{clade:11} {d['d_nodes']:>7} {dm:>6} {status:7} {d['relation']}")
        if d["diff_muts"]:
            lines.append(f"{'':11} mutations between anchors: {d['diff_muts']}")
    lines.append("")
    lines.append(f"max observed patristic distance: {worst} mutation(s) "
                 f"(threshold {max_dist}).")
    lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_distance_data(path, distances):
    """Machine-readable distances for the check_anchors rule."""
    doc = {}
    for clade, d in distances.items():
        doc[clade] = {
            "rep_absent": bool(d.get("rep_absent")),
            "error": d.get("error"),
            "d_nodes": d.get("d_nodes"),
            "d_muts": d.get("d_muts"),
            "relation": d.get("relation"),
            "diff_muts": d.get("diff_muts", []),
        }
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=True, default_flow_style=False)


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml")
    p.add_argument("--nodes-representative")
    p.add_argument("--nodes-founder")
    p.add_argument("--sample-paths")
    p.add_argument("--clade-paths")
    p.add_argument("--out-yaml", help="Nextstrain-clade output YAML.")
    p.add_argument("--out-pango", help="PANGO-clade output YAML.")
    p.add_argument("--out-comparison", help="Human anchor-distance report.")
    p.add_argument("--out-distances", help="Machine-readable anchor distances.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    log_file = cfg["log"]
    lib.validate_resolution_overrides(cfg["representative_name_overrides"],
                                      cfg["representative_absent_ok"])

    nodes_rep = cli.nodes_representative or cfg["nodes_representative"]
    nodes_founder = cli.nodes_founder or cfg["nodes_founder"]
    sample_paths = cli.sample_paths
    clade_paths = cli.clade_paths
    out_yaml = cli.out_yaml or cfg["output"]
    out_pango = cli.out_pango or cfg["output_pango"]
    out_comparison = cli.out_comparison or cfg["comparison_report"]
    out_distances = cli.out_distances or cfg["anchor_distances"]
    if not sample_paths:
        raise ValueError("need --sample-paths")
    if not clade_paths:
        raise ValueError("need --clade-paths")

    for path in (out_yaml, out_pango, out_comparison, out_distances):
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

    # The output is founder-anchored: every clade must have a founder node.
    founder_node, _, no_founder = lib.read_nodes(nodes_founder)
    if no_founder:
        raise ValueError(
            f"{len(no_founder)} PANGO clade(s) have no annotated founder node and "
            f"cannot be anchored: {sorted(no_founder)}")
    rep_node, _, rep_unresolved = lib.read_nodes(nodes_rep)
    lib.log(f"Founder anchors: {len(founder_node)}; representative anchors: "
            f"{len(rep_node)} ({len(rep_unresolved)} absent)", log_file)

    # --- Anchor cross-check table (report-only here; enforced in check_anchors) --
    lib.log("\n=== Cross-checking founder vs representative anchors ===", log_file)
    distances = lib.compute_anchor_distances(founder_node, rep_node, sample_paths,
                                             log_file)
    write_distance_report(out_comparison, distances, cfg["max_anchor_distance"],
                          set(cfg["representative_absent_ok"]))
    write_distance_data(out_distances, distances)
    lib.log(f"Wrote anchor distance table to {out_comparison} / {out_distances}",
            log_file)

    # --- Build both outputs from the founder anchoring (one shared scan) --------
    ns_founder_nodes = {clade_founders[v] for v in ns_full_values
                        if v in clade_founders}
    shared_targets = set(founder_node.values()) | ns_founder_nodes
    lib.log("\n=== Capturing ancestor chains (shared across both outputs) ===",
            log_file)
    chains = lib.capture_ancestor_chains(sample_paths, shared_targets, log_file)

    lib.log("\n=== Building Nextstrain-clade masks (founder-anchored) ===", log_file)
    ns_masks = lib.build_masks_for_nodes(
        spec, founder_node, ns_full_values, clade_founders, sample_paths,
        start, end, log_file, chains=chains)

    lib.log("\n=== Building PANGO-clade masks (founder-anchored) ===", log_file)
    pango_masks = lib.build_masks_for_nodes(
        spec, founder_node, list(spec), clade_founders, sample_paths,
        start, end, log_file, chains=chains)

    # The output is written as-is; the check_recovery rule validates old-scheme
    # recovery and never edits the output to fill a gap.
    dropped_exclusions = {c: e["exclusions"] for c, e in spec.items() if e["exclusions"]}
    tree_name = os.path.basename(cfg["tree"]) if cfg["tree"] else "(sample-paths supplied)"
    base_prov = [
        f"source spec:   {cfg['branch_spec']}",
        f"tree:          {tree_name}",
        f"metadata:      {os.path.basename(cfg['metadata'])}",
        f"clade column:  {cfg['clade_column']}",
        f"include range: [{start}, {end}] (masks outside this range are omitted)",
        "anchoring:     founder (MRCA) node per PANGO clade",
    ]
    lib.write_output(out_yaml, ns_masks, spec, base_prov, dropped_exclusions,
                     level="Nextstrain")
    lib.log(f"Wrote {out_yaml}", log_file)
    lib.write_output(out_pango, pango_masks, spec, base_prov, dropped_exclusions,
                     level="PANGO")
    lib.log(f"Wrote {out_pango}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
