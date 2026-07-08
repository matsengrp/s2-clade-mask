#!/usr/bin/env python3
"""Generate the Nextstrain-clade-keyed ``usher_masked_sites.yaml`` in one shot.

This is a thin convenience wrapper around the :mod:`s2clademask` library that runs
every stage in a single process -- handy for a quick local run. The staged
Snakemake pipeline (``Snakefile``) is the primary path: it runs the exact same
library functions as separate, cached rules (the two matUtils passes, the two
anchoring-mode resolutions, the comparison + build, and the old-scheme recovery
check), so troubleshooting does not re-pay for the ~1 h / ~19 GB ``--sample-paths``
extract every cycle.

Both paths produce identical results because both call the same code:
representative-primary resolution with human-approved per-clade fallbacks
(design invariant #1), exact name matching only (#2), config-approved deviations
(#3), superset recovery of the vendored old scheme (#4), and UTR drop (#5). See
``README.md`` for the full algorithm write-up.

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
                                      args.representative_fallbacks)

    # Ensure the scratch dir (for the large paths files) and output dir exist.
    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    lib.log(f"Parsing branch spec: {args.branch_spec}", log_file)
    spec, version = lib.parse_branch_spec(args.branch_spec)
    lib.log(f"  {len(spec)} PANGO clades, spec version {version}", log_file)

    # Obtain sample paths (root->tip node IDs for every sample; used for the
    # membership scan and as the representative-sample resolver).
    tmp_paths = None
    sample_paths = args.sample_paths
    if sample_paths is None:
        if not args.tree:
            raise ValueError("config needs 'tree' (or a precomputed 'sample_paths')")
        fd, tmp_paths = tempfile.mkstemp(prefix="sample_paths_", suffix=".tsv",
                                         dir=args.workdir)
        os.close(fd)
        sample_paths = lib.run_sample_paths(args.tree, tmp_paths, args.matutils, log_file)

    # Obtain clade paths (annotated founder node per clade).
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

    try:
        lib.log("Resolving PANGO clade nodes (representative sample primary; "
                "approved overrides/fallbacks only)...", log_file)
        clade_to_node, _via = lib.resolve_clade_nodes(
            spec, clade_founders, sample_paths,
            args.representative_name_overrides, args.representative_fallbacks,
            log_file)

        lib.log(f"Loading Nextstrain clades from {args.metadata}...", log_file)
        ns_full_values = lib.load_nextstrain_clade_values(
            args.metadata, args.clade_column, log_file)

        lib.log("Assigning PANGO masks to Nextstrain clades "
                "(ancestor-or-nested over-mask; PANGO clades anchored at their "
                "representative node)...", log_file)
        clade_masks = lib.build_masks_for_nodes(
            spec, clade_to_node, ns_full_values, clade_founders, sample_paths,
            args.site_include_start, args.site_include_end, log_file)
    finally:
        if tmp_paths and os.path.exists(tmp_paths):
            os.remove(tmp_paths)
        if tmp_clade_paths and os.path.exists(tmp_clade_paths):
            os.remove(tmp_clade_paths)

    # Report which PANGO clades contribute masks to each Nextstrain clade. When
    # a Nextstrain clade lists PANGO clades finer than itself, its mask has been
    # widened to the whole clade (over-mask coarsening).
    lib.log("Mask sources per Nextstrain clade (contributing PANGO clades):", log_file)
    for ns_clade in sorted(clade_masks):
        sources = [s for s in clade_masks[ns_clade]["sources"]
                   if not s.startswith("old:")]
        if sources:
            lib.log(f"  {ns_clade} <- {sources}", log_file)

    dropped_exclusions = {c: e["exclusions"] for c, e in spec.items() if e["exclusions"]}

    missing = {}
    if args.old_yaml:
        old_sites = lib.load_old_sites(args.old_yaml)
        missing = lib.check_recovers_old(
            old_sites, clade_masks, args.union_old,
            include_range=(args.site_include_start, args.site_include_end),
            log_file=log_file)

    provenance = [
        f"source spec:   {args.branch_spec}",
        f"tree:          {os.path.basename(args.tree) if args.tree else '(sample-paths supplied)'}",
        f"metadata:      {os.path.basename(args.metadata)}",
        f"clade column:  {args.clade_column}",
        f"include range: [{args.site_include_start}, {args.site_include_end}] "
        "(masks outside this range are omitted)",
    ]
    lib.write_output(args.output, clade_masks, spec, provenance, dropped_exclusions)
    lib.log(f"Wrote {args.output}", log_file)

    if missing and args.strict and not args.union_old:
        lib.log("ERROR: --strict and old masks not fully recovered (see above).",
                log_file)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
