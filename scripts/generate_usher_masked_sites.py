#!/usr/bin/env python3
"""Generate the founder-anchored masking YAMLs in one shot.

Thin convenience wrapper around the :mod:`s2clademask` library that runs every
stage in a single process. The staged Snakemake pipeline (``Snakefile``) is the
primary path: it runs the exact same library functions as separate, cached rules
(the two matUtils passes, the two anchoring-mode resolutions, the build + anchor
cross-check, and the old-scheme recovery check), so troubleshooting does not
re-pay for the ~1 h / ~19 GB ``--sample-paths`` extract every cycle.

Both paths produce identical results because both call the same code. Each PANGO
clade is anchored at its annotated **founder (MRCA) node** (design invariant #1);
the **representative sample** is used only to cross-check that the founder sits
within ``max_anchor_distance`` mutations of it (else a hard error, unless the
representative is an approved-absent one). Two YAMLs are written: the
Nextstrain-clade ``output`` (for mat-pcp) and the per-PANGO-lineage
``output_pango`` companion. See ``README.md`` for the full write-up.

Usage (from the repo root so relative config paths resolve):

    python scripts/generate_usher_masked_sites.py --config config.yaml
"""

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml",
                   help="Path to the YAML config (default: config.yaml). See the "
                        "repo's config.yaml for the schema.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    args = argparse.Namespace(**cfg)
    log_file = args.log
    lib.log(f"Loaded config: {cli.config}", log_file)
    lib.validate_resolution_overrides(args.representative_name_overrides,
                                      args.representative_absent_ok)

    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)
    for out in (args.output, args.output_pango):
        out_dir = os.path.dirname(out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    lib.log(f"Parsing branch spec: {args.branch_spec}", log_file)
    spec, version = lib.parse_branch_spec(args.branch_spec)
    lib.log(f"  {len(spec)} PANGO clades, spec version {version}", log_file)

    # Sample paths (root->tip node IDs; membership scan + representative resolver).
    tmp_paths = None
    sample_paths = args.sample_paths
    if sample_paths is None:
        if not args.tree:
            raise ValueError("config needs 'tree' (or a precomputed 'sample_paths')")
        fd, tmp_paths = tempfile.mkstemp(prefix="sample_paths_", suffix=".tsv",
                                         dir=args.workdir)
        os.close(fd)
        sample_paths = lib.run_sample_paths(args.tree, tmp_paths, args.matutils, log_file)

    # Clade paths (annotated founder node per clade; the primary anchor).
    tmp_clade_paths = None
    clade_paths = args.clade_paths
    if clade_paths is None:
        if not args.tree:
            raise ValueError("config needs 'tree' (or a precomputed 'clade_paths')")
        fd, tmp_clade_paths = tempfile.mkstemp(prefix="clade_paths_", suffix=".tsv",
                                               dir=args.workdir)
        os.close(fd)
        clade_paths = lib.run_clade_paths(args.tree, tmp_clade_paths, args.matutils,
                                          log_file)
    clade_founders = lib.parse_clade_paths(clade_paths)
    lib.log(f"Parsed {len(clade_founders):,} clade annotations from {clade_paths}",
            log_file)

    anchor_failures = []
    try:
        lib.log("Resolving founder (MRCA) anchor for every PANGO clade...", log_file)
        founder_node = lib.resolve_founder_or_fail(spec, clade_founders, log_file)

        lib.log("Resolving representative anchors (for the founder cross-check)...",
                log_file)
        rep_node, _via, rep_unresolved = lib.resolve_representative_only(
            spec, sample_paths, args.representative_name_overrides, log_file)
        not_ok = [c for c in rep_unresolved
                  if c not in args.representative_absent_ok]
        if not_ok:
            raise ValueError(
                "Representative absent from the tree and not in "
                f"representative_absent_ok for clade(s): {sorted(not_ok)}")

        lib.log("Cross-checking founder vs representative anchors...", log_file)
        distances = lib.compute_anchor_distances(founder_node, rep_node,
                                                 sample_paths, log_file)
        for clade in sorted(distances):
            d = distances[clade]
            if d.get("rep_absent"):
                continue
            if d.get("error"):
                anchor_failures.append(f"{clade}: {d['error']}")
            elif d["d_muts"] > args.max_anchor_distance:
                anchor_failures.append(
                    f"{clade}: founder {d['d_muts']} muts from representative "
                    f"(> {args.max_anchor_distance}); {d['diff_muts']}")
            else:
                lib.log(f"  {clade}: {d['d_muts']} mutation(s) between anchors "
                        f"({d['relation']})", log_file)

        lib.log(f"Loading Nextstrain clades from {args.metadata}...", log_file)
        ns_full_values = lib.load_nextstrain_clade_values(
            args.metadata, args.clade_column, log_file)

        start, end = args.site_include_start, args.site_include_end
        ns_founder_nodes = {clade_founders[v] for v in ns_full_values
                            if v in clade_founders}
        chains = lib.capture_ancestor_chains(
            sample_paths, set(founder_node.values()) | ns_founder_nodes, log_file)

        lib.log("Building Nextstrain-clade masks (founder-anchored)...", log_file)
        ns_masks = lib.build_masks_for_nodes(
            spec, founder_node, ns_full_values, clade_founders, sample_paths,
            start, end, log_file, chains=chains)
        lib.log("Building PANGO-clade masks (founder-anchored)...", log_file)
        pango_masks = lib.build_masks_for_nodes(
            spec, founder_node, list(spec), clade_founders, sample_paths,
            start, end, log_file, chains=chains)
    finally:
        if tmp_paths and os.path.exists(tmp_paths):
            os.remove(tmp_paths)
        if tmp_clade_paths and os.path.exists(tmp_clade_paths):
            os.remove(tmp_clade_paths)

    missing = {}
    if args.old_yaml:
        old_sites = lib.load_old_sites(args.old_yaml)
        rng = (args.site_include_start, args.site_include_end)
        missing = lib.check_recovers_old(
            old_sites, ns_masks, args.union_old, include_range=rng,
            log_file=log_file)
        # Per-clade added-mask summary table (regenerated each run; mirrors README).
        lib.write_recovery_summary(args.recovery_summary, old_sites, ns_masks, rng,
                                   gen_yaml=args.output, old_yaml=args.old_yaml)
        lib.log(f"Wrote added-mask summary to {args.recovery_summary}", log_file)

    dropped_exclusions = {c: e["exclusions"] for c, e in spec.items() if e["exclusions"]}
    tree_name = os.path.basename(args.tree) if args.tree else "(sample-paths supplied)"
    provenance = [
        f"source spec:   {args.branch_spec}",
        f"tree:          {tree_name}",
        f"metadata:      {os.path.basename(args.metadata)}",
        f"clade column:  {args.clade_column}",
        f"include range: [{args.site_include_start}, {args.site_include_end}] "
        "(masks outside this range are omitted)",
        "anchoring:     founder (MRCA) node per PANGO clade",
    ]
    lib.write_output(args.output, ns_masks, spec, provenance, dropped_exclusions,
                     level="Nextstrain")
    lib.log(f"Wrote {args.output}", log_file)
    lib.write_output(args.output_pango, pango_masks, spec, provenance,
                     dropped_exclusions, level="PANGO")
    lib.log(f"Wrote {args.output_pango}", log_file)

    rc = 0
    if anchor_failures:
        lib.log("ERROR: anchor cross-check failed:", log_file)
        for msg in anchor_failures:
            lib.log(f"  {msg}", log_file)
        rc = 1
    if missing and args.strict and not args.union_old:
        lib.log("ERROR: --strict and old masks not fully recovered (see above).",
                log_file)
        rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
