#!/usr/bin/env python3
"""Pipeline rule entrypoint: check the generated YAML recovers the old scheme.

Reads the written ``usher_masked_sites.yaml`` back (``load_generated_masks``) and
verifies it is a **superset** of the vendored old hand-made scheme within the
``site_include_range`` (``check_recovers_old``): every old ``(clade, site)`` mask
must still be masked as a site, or now at the mutation level (an intentional
site->mutation downgrade). Old UTR masks outside the range are reported as
intentionally-dropped, not failures.

Writes the pass/fail report to ``--out`` and, alongside it, a per-clade
added-mask summary table (``--out-summary`` / config ``recovery_summary``) that
quantifies how many masked positions the new scheme adds over the old one --
regenerated every run so the README snapshot can be refreshed. Under
``strict: true`` in the config, exits nonzero on a *true* shortfall (design
invariant #4 -- a real shortfall is a bug to investigate, not to accept).

The recovery logic lives in this script (it is used by no other rule); config
loading and the IO helpers come from ``s2clademask``.

    python scripts/check_recovery.py -c config.yaml \\
        --yaml results/usher_masked_sites.yaml --out results/recovery_report.txt \\
        --out-summary results/recovery_summary.txt
"""

import argparse
import os
import sys
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


# ---------------------------------------------------------------------------
# Load the old scheme and the freshly-generated output
# ---------------------------------------------------------------------------


def load_old_sites(old_yaml):
    """Return ``{clade: set(sites)}`` from an existing usher_masked_sites.yaml."""
    with lib.open_maybe_gzip(old_yaml) as f:
        old = yaml.safe_load(f)
    # Blocks are keyed by their ``clades`` list; one block can name several clades.
    result = defaultdict(set)
    for _, block in old.items():
        clades = block.get("clades", []) or []
        sites = set(block.get("sites", []) or [])
        for c in clades:
            result[c].update(sites)
    return result


def load_generated_masks(yaml_path):
    """Read a written ``usher_masked_sites.yaml`` into the in-memory mask map.

    Returns ``{bare_clade: {"sites": set[int], "mutations": set[str]}}`` -- the
    same shape ``build_clade_masks`` produces (minus ``sources``) -- so this
    standalone rule can re-check the on-disk output without re-running the
    assignment. Blocks are keyed by their ``clades`` list, matching how
    ``load_old_sites`` reads the old scheme.
    """
    with lib.open_maybe_gzip(yaml_path) as f:
        doc = yaml.safe_load(f) or {}
    result = {}
    for _, block in doc.items():
        clades = block.get("clades", []) or []
        sites = set(int(s) for s in (block.get("sites", []) or []))
        muts = set(block.get("mutations", []) or [])
        for c in clades:
            entry = result.setdefault(c, {"sites": set(), "mutations": set()})
            entry["sites"].update(sites)
            entry["mutations"].update(muts)
    return result


# ---------------------------------------------------------------------------
# Superset recovery check
# ---------------------------------------------------------------------------


def check_recovers_old(old_sites, clade_masks, include_range=None, log_file=None):
    """Verify the new file recovers every old ``(clade, site)`` mask.

    An old site is recovered if it is still masked as a **site**, or if it is now
    masked at the **mutation** level (>=1 new mutation at that position for the
    clade). The latter is an intentional, more-precise "downgrade" (the old file
    masked whole sites; the new UCSC spec masks specific reversion mutations).
    Only sites recovered by *neither* mechanism are true failures.

    If ``include_range`` is given as ``(start, end)``, old sites outside it (the
    UTRs) are **intentionally not recovered** -- they are dropped from the new
    scheme because the pipeline's range filter discards them anyway -- and are
    reported separately rather than counted as failures.

    Returns ``truly_missing`` (``{clade: set(sites)}``).
    """
    truly_missing = defaultdict(set)
    downgraded = defaultdict(set)
    out_of_range = defaultdict(set)
    lo, hi = include_range if include_range else (None, None)
    for clade, sites in old_sites.items():
        block = clade_masks.get(clade, {})
        new_site_set = block.get("sites", set())
        # Positions now covered at the mutation level count as recovered, too.
        new_mut_positions = {int(m[1:-1]) for m in block.get("mutations", set())}
        for site in sites:
            if include_range and not (lo <= site <= hi):
                out_of_range[clade].add(site)
                continue
            if site in new_site_set:
                continue
            if site in new_mut_positions:
                downgraded[clade].add(site)
            else:
                truly_missing[clade].add(site)

    if out_of_range:
        n = sum(len(v) for v in out_of_range.values())
        lib.log(f"OLD-SCHEME RECOVERY: {n} old (clade, site) mask(s) fall outside the "
                f"site_include_range [{lo}, {hi}] and are intentionally NOT recovered "
                "(UTR sites dropped; discarded by the range filter at runtime).",
                log_file)

    if downgraded:
        n = sum(len(v) for v in downgraded.values())
        lib.log(f"OLD-SCHEME RECOVERY: {n} old whole-site mask(s) are now masked at "
                "the mutation level (intentional site->mutation downgrade):", log_file)
        for clade in sorted(downgraded):
            lib.log(f"  {clade}: {sorted(downgraded[clade])}", log_file)

    if truly_missing:
        lib.log("OLD-SCHEME RECOVERY: the following (clade, site) masks are present "
                "in the old file but NOT recovered (neither site nor mutation):",
                log_file)
        for clade in sorted(truly_missing):
            lib.log(f"  {clade}: {sorted(truly_missing[clade])}", log_file)
    else:
        lib.log("OLD-SCHEME RECOVERY: OK -- every old (clade, site) mask is recovered "
                "as a site or mutation mask.", log_file)
    return truly_missing


# ---------------------------------------------------------------------------
# Per-clade added-mask summary (informational; regenerated every run)
# ---------------------------------------------------------------------------


def summarize_added_masks(old_sites, clade_masks, include_range):
    """Per-clade tally of masked positions the new scheme adds over the old one.

    A genome position counts as *masked* for a clade whether it is masked
    whole-site or via a mutation. Restricted to ``include_range`` (old UTR masks are
    dropped by design, so counting them would be misleading). Returns
    ``(rows, totals)`` where ``rows`` is a list of dicts sorted by clade
    (``{clade, old, new, added, new_clade}``) and ``totals`` carries the aggregate
    counts (grand totals, split of added masks between old-scheme-absent clades and
    the shared clades, and distinct-position growth across all clades).
    """
    lo, hi = include_range

    def positions(sites, muts=()):
        # Collapse whole-site and mutation-level masks to the set of masked positions.
        pos = {int(s) for s in sites if lo <= int(s) <= hi}
        for m in muts:
            p = int(m[1:-1])
            if lo <= p <= hi:
                pos.add(p)
        return pos

    old_pos = {c: positions(s) for c, s in old_sites.items()}
    new_pos = {c: positions(b.get("sites", ()), b.get("mutations", ()))
               for c, b in clade_masks.items()}
    old_clades = set(old_pos)

    rows = []
    added_in_new_clades = added_in_shared = 0
    for clade in sorted(set(old_pos) | set(new_pos)):
        o, n = old_pos.get(clade, set()), new_pos.get(clade, set())
        added = n - o
        is_new = clade not in old_clades  # clade absent from the old hand-made scheme
        rows.append({"clade": clade, "old": len(o), "new": len(n),
                     "added": len(added), "new_clade": is_new})
        if is_new:
            added_in_new_clades += len(added)
        else:
            added_in_shared += len(added)

    old_all = set().union(*old_pos.values()) if old_pos else set()
    new_all = set().union(*new_pos.values()) if new_pos else set()
    totals = {
        "old": sum(r["old"] for r in rows),
        "new": sum(r["new"] for r in rows),
        "added": sum(r["added"] for r in rows),
        "added_in_new_clades": added_in_new_clades,
        "added_in_shared": added_in_shared,
        "n_new_clades": sum(1 for r in rows if r["new_clade"]),
        "distinct_old": len(old_all),
        "distinct_new": len(new_all),
        "distinct_added": len(new_all - old_all),
    }
    return rows, totals


def write_recovery_summary(out_path, old_sites, clade_masks, include_range,
                           gen_yaml=None, old_yaml=None):
    """Write the per-clade added-masks summary table (regenerated each run).

    A companion to the pass/fail ``recovery_report``: it quantifies, per Nextstrain
    clade, how many in-range masked positions the new scheme carries vs the old
    hand-made scheme and how many it adds. Mirrors the snapshot table in README.md
    so that snapshot can be refreshed from a real run.
    """
    lo, hi = include_range
    rows, t = summarize_added_masks(old_sites, clade_masks, include_range)
    lines = []
    lines.append("=" * 66)
    lines.append("ADDED-MASK SUMMARY  (new scheme  vs  old hand-made scheme)")
    lines.append(f"masked positions in range [{lo}, {hi}]; a position counts whether")
    lines.append("masked whole-site or via a mutation. 'added' = masked now, not before.")
    if gen_yaml or old_yaml:
        lines.append(f"new: {gen_yaml}")
        lines.append(f"old: {old_yaml}")
    lines.append("=" * 66)
    lines.append("")
    lines.append(f"{'clade':9} {'old':>5} {'new':>5} {'added':>6}  note")
    for r in rows:
        note = "new clade (absent from old scheme)" if r["new_clade"] else ""
        lines.append(f"{r['clade']:9} {r['old']:>5} {r['new']:>5} "
                     f"{r['added']:>6}  {note}")
    lines.append("-" * 66)
    lines.append(f"{'TOTAL':9} {t['old']:>5} {t['new']:>5} {t['added']:>6}")
    lines.append("")
    lines.append(f"Added masks: {t['added']} total "
                 f"({t['added_in_new_clades']} in the {t['n_new_clades']} clades "
                 f"absent from the old scheme, {t['added_in_shared']} in the shared "
                 "clades the old scheme already had).")
    lines.append(f"Distinct genomic positions masked (any clade): "
                 f"{t['distinct_old']} -> {t['distinct_new']} "
                 f"(+{t['distinct_added']} never masked in any clade before).")
    lines.append("")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return rows, t


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml")
    p.add_argument("--yaml", help="Generated usher_masked_sites.yaml to check.")
    p.add_argument("--old", help="Old-scheme baseline YAML (default: config old_yaml).")
    p.add_argument("--out", help="Recovery report text file.")
    p.add_argument("--out-summary", help="Per-clade added-mask summary table.")
    cli = p.parse_args(argv)

    cfg = lib.load_config(cli.config)
    gen_yaml = cli.yaml or cfg["output"]
    old_yaml = cli.old or cfg["old_yaml"]
    out_report = cli.out or cfg["recovery_report"]
    out_summary = cli.out_summary or cfg["recovery_summary"]
    start, end = cfg["site_include_start"], cfg["site_include_end"]

    for path in (out_report, out_summary):
        d = os.path.dirname(os.path.abspath(path))
        if d:
            os.makedirs(d, exist_ok=True)

    # Log to the report file (check_recovers_old logs its findings line by line).
    open(out_report, "w").close()

    if not old_yaml:
        lib.log("No 'old_yaml' configured -- skipping old-scheme recovery check.",
                out_report)
        with open(out_summary, "w") as f:
            f.write("No 'old_yaml' configured -- added-mask summary not computed.\n")
        return 0

    lib.log(f"OLD-SCHEME RECOVERY CHECK: {gen_yaml} vs {old_yaml} "
            f"(range [{start}, {end}], strict={cfg['strict']})", out_report)
    clade_masks = load_generated_masks(gen_yaml)
    old_sites = load_old_sites(old_yaml)
    # Validate only -- the on-disk output is never mutated to paper over a gap.
    truly_missing = check_recovers_old(
        old_sites, clade_masks, include_range=(start, end), log_file=out_report)

    # Per-clade added-mask summary table (regenerated each run; mirrors README).
    _, totals = write_recovery_summary(
        out_summary, old_sites, clade_masks, (start, end),
        gen_yaml=gen_yaml, old_yaml=old_yaml)
    lib.log(f"OLD-SCHEME RECOVERY: wrote added-mask summary ({totals['added']} "
            f"masks added over the old scheme) to {out_summary}", out_report)

    if truly_missing and cfg["strict"]:
        lib.log("ERROR: strict=true and old masks not fully recovered (see above).",
                out_report)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
