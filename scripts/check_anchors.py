#!/usr/bin/env python3
"""Pipeline rule entrypoint: enforce the founder-vs-representative closeness assertion.

Reads the machine-readable anchor-distance table written by ``compare_and_build``
and fails (nonzero exit) if any clade's founder anchor sits more than
``max_anchor_distance`` mutations (patristic) from its representative anchor, or if
a clade's representative is absent from the tree without being listed in
``representative_absent_ok``. This is a hard guardrail: it keeps the deliberate
divergence from upstream's representative-node choice provably tiny.

    python scripts/check_anchors.py -c config.yaml \\
        --distances results/anchor_distances.yaml --out results/anchor_check.txt
"""

import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml")
    p.add_argument("--distances", help="anchor_distances.yaml from compare_and_build.")
    p.add_argument("--out", help="Assertion-result report.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    distances_path = cli.distances or cfg["anchor_distances"]
    out_report = cli.out or cfg["anchor_check_report"]
    max_dist = cfg["max_anchor_distance"]
    absent_ok = set(cfg["representative_absent_ok"])

    out_dir = os.path.dirname(os.path.abspath(out_report))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    open(out_report, "w").close()  # truncate; lib.log appends line by line below

    # The per-clade distance table produced by compare_and_build.
    with lib.open_maybe_gzip(distances_path) as f:
        distances = yaml.safe_load(f) or {}

    lib.log(f"ANCHOR CHECK: max_anchor_distance = {max_dist} mutations (patristic); "
            f"{len(distances)} clade(s)", out_report)

    # Collect every violation and report them together, rather than failing on the
    # first, so one run surfaces all anchors that need attention.
    failures = []
    for clade in sorted(distances):
        d = distances[clade]
        if d.get("rep_absent"):
            # An absent representative is only OK if the maintainer pre-approved it.
            if clade not in absent_ok:
                failures.append(f"{clade}: representative ABSENT from tree and not in "
                                "representative_absent_ok")
            else:
                lib.log(f"  {clade}: representative absent (approved) -- skipped",
                        out_report)
            continue
        if d.get("error"):
            failures.append(f"{clade}: {d['error']}")
            continue
        dm = d.get("d_muts")
        if dm is None:
            failures.append(f"{clade}: missing distance data")
        elif dm > max_dist:
            failures.append(f"{clade}: founder is {dm} mutations from representative "
                            f"(> {max_dist}); mutations between: {d.get('diff_muts')}")
        else:
            lib.log(f"  {clade}: {dm} mutation(s) -- OK", out_report)

    if failures:
        lib.log("\nANCHOR CHECK FAILED:", out_report)
        for msg in failures:
            lib.log(f"  {msg}", out_report)
        lib.log("\nInvestigate by hand: either the founder/representative genuinely "
                "diverged (a spec or tree change), or add the clade to "
                "representative_absent_ok / raise max_anchor_distance only after "
                "confirming the anchors really are close.", out_report)
        return 1

    lib.log("\nANCHOR CHECK OK -- every founder anchor is within "
            f"{max_dist} mutation(s) of its representative (or approved absent).",
            out_report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
