#!/usr/bin/env python3
"""Pipeline rule entrypoint: check the generated YAML recovers the old scheme.

Reads the written ``usher_masked_sites.yaml`` back (``load_generated_masks``) and
verifies it is a **superset** of the vendored old hand-made scheme within the
``site_include_range`` (``check_recovers_old``): every old ``(clade, site)`` mask
must still be masked as a site, or now at the mutation level (an intentional
site->mutation downgrade). Old UTR masks outside the range are reported as
intentionally-dropped, not failures.

Writes the report to ``--out`` and, under ``strict: true`` in the config, exits
nonzero on a *true* shortfall (design invariant #4 -- a real shortfall is a bug to
investigate, not to accept).

    python scripts/check_recovery.py -c config.yaml \\
        --yaml results/usher_masked_sites.yaml --out results/recovery_report.txt
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml")
    p.add_argument("--yaml", help="Generated usher_masked_sites.yaml to check.")
    p.add_argument("--old", help="Old-scheme baseline YAML (default: config old_yaml).")
    p.add_argument("--out", help="Recovery report text file.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    gen_yaml = cli.yaml or cfg["output"]
    old_yaml = cli.old or cfg["old_yaml"]
    out_report = cli.out or cfg["recovery_report"]
    start, end = cfg["site_include_start"], cfg["site_include_end"]

    out_dir = os.path.dirname(os.path.abspath(out_report))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Log to the report file (check_recovers_old logs its findings line by line).
    open(out_report, "w").close()

    if not old_yaml:
        lib.log("No 'old_yaml' configured -- skipping old-scheme recovery check.",
                out_report)
        return 0

    lib.log(f"OLD-SCHEME RECOVERY CHECK: {gen_yaml} vs {old_yaml} "
            f"(range [{start}, {end}], strict={cfg['strict']})", out_report)
    clade_masks = lib.load_generated_masks(gen_yaml)
    old_sites = lib.load_old_sites(old_yaml)
    # union_old is applied at build time (compare_and_build); here we only report
    # and, under strict, gate -- so never mutate the on-disk output.
    truly_missing = lib.check_recovers_old(
        old_sites, clade_masks, union_old=False, include_range=(start, end),
        log_file=out_report)

    if truly_missing and cfg["strict"] and not cfg["union_old"]:
        lib.log("ERROR: strict=true and old masks not fully recovered (see above).",
                out_report)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
