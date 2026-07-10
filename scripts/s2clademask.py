#!/usr/bin/env python3
"""Shared helpers for the s2-clade-mask generator.

This module holds only the pieces the staged Snakemake pipeline shares across
**more than one** rule script: small IO utilities (``open_maybe_gzip``, ``log``),
the two input parsers (``parse_branch_spec`` for the UCSC spec, ``parse_clade_paths``
for the ``--clade-paths`` TSV), the config-supplied override validation
(``validate_resolution_overrides``), the node-map artifact schema exchanged between
rules (``write_nodes`` / ``read_nodes``), and config loading (``load_config``).

Logic used by a **single** rule lives in that rule's script, not here:

* ``scripts/resolve_nodes.py``     -- representative / founder node resolution
* ``scripts/compare_and_build.py`` -- the over-mask assignment, anchor cross-check,
                                      and YAML emission
* ``scripts/check_recovery.py``    -- old-scheme recovery + added-mask summary

The two ``matUtils extract`` passes are run directly by the Snakefile's extract
rules. The full algorithm write-up lives in ``README.md``; the design invariants
(founder-anchored resolution with a representative cross-check, exact matching only,
human-approved deviations, superset recovery, UTR drop) are documented there and in
``CLAUDE.md``.
"""

import re
import gzip

import yaml


def open_maybe_gzip(path, mode="rt"):
    """Open ``path`` transparently whether or not it is gzip-compressed."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def log(msg, log_file=None):
    """Print ``msg`` and, if given, append it to ``log_file``."""
    print(msg, flush=True)
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")


# ---------------------------------------------------------------------------
# 1. Parse the UCSC spec (shared: resolve_nodes + compare_and_build)
# ---------------------------------------------------------------------------

# Reversion mutations are directional refPOSalt strings, e.g. "G670T".
REVERSION_RE = re.compile(r"^[ACGT][0-9]+[ACGT]$")


def parse_branch_spec(path):
    """Parse ``branchSpecificMask.yml``.

    Returns ``(spec, version)`` where ``spec`` maps each PANGO clade name to a
    dict with keys: ``representative`` (str), ``backtrack`` (int),
    ``sites`` (set[int], ranges expanded and merged with individual sites),
    ``mutations`` (list[str] reversions), ``exclusions`` (list[dict]).
    """
    with open_maybe_gzip(path) as f:
        raw = yaml.safe_load(f)

    # The spec carries a top-level ``version`` key alongside the per-clade blocks;
    # pop it out so the remaining keys are all clades. A missing version means the
    # file is not the spec we expect, so fail loudly rather than silently.
    version = raw.pop("version", None)
    if version is None:
        raise ValueError(f"No 'version' key found in {path}")

    spec = {}
    for clade, entry in raw.items():
        rep = entry.get("representative")
        if not rep:
            raise ValueError(f"No 'representative' for clade {clade} in {path}")

        # ``ranges: [[lo, hi]]`` are inclusive integer spans; expand and merge them
        # with any individual ``sites`` into one flat site set.
        sites = set()
        for lo, hi in entry.get("ranges", []) or []:
            sites.update(range(int(lo), int(hi) + 1))
        for s in entry.get("sites", []) or []:
            sites.add(int(s))

        # Validate each reversion string up front so a malformed entry fails here,
        # not deep in the mask-assignment stage.
        mutations = []
        for rev in entry.get("reversions", []) or []:
            if not REVERSION_RE.match(rev):
                raise ValueError(f"Malformed reversion '{rev}' for clade {clade}")
            mutations.append(rev)

        if not sites and not mutations:
            # A clade with only exclusions and nothing to mask is unusual; keep
            # it so its representative still anchors ancestry, but warn.
            pass

        spec[clade] = {
            "representative": rep,
            "backtrack": int(entry.get("representativeBacktrack", 0) or 0),
            "sites": sites,
            "mutations": mutations,
            "exclusions": entry.get("exclusions", []) or [],
        }
    return spec, version


# ---------------------------------------------------------------------------
# 2. Parse the clade-paths TSV (shared: resolve_nodes + compare_and_build)
#
# The two ``matUtils extract`` passes (--sample-paths / --clade-paths) are run
# directly by the Snakefile's extract rules; see the note there about the ``-d /``
# absolute-output-path quirk. This module only parses the TSVs they produce.
# ---------------------------------------------------------------------------


def parse_clade_paths(path):
    """Parse ``matUtils --clade-paths`` output into ``{clade: founder_node_id}``.

    The file is TSV with header ``clade<TAB>root_id<TAB>from_tree_root``; the
    ``root_id`` column is the node ID of the clade's founder (MRCA). The tree is
    annotated with both Nextstrain clades and PANGO lineages, so PANGO clade
    names from the spec (e.g. ``BA.4``, ``BA.2.75``) resolve directly.
    """
    founders = {}
    with open_maybe_gzip(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # setdefault: keep the first row for a clade if it appears more than once.
            if len(parts) >= 2 and parts[0] and parts[1]:
                founders.setdefault(parts[0], parts[1])
    return founders


# ---------------------------------------------------------------------------
# Human-approved resolution overrides (shared: resolve_nodes + compare_and_build)
#
# Each PANGO clade is anchored at its annotated **founder (MRCA) node** -- always
# present in a full gisaid-and-public MAT, and a hard error if missing. The clade's
# **representative sample** is resolved only to *cross-check* the founder: we assert
# the representative and founder anchors are within ``max_anchor_distance``
# mutations of each other (patristic). Representative samples are matched to the
# tree by EXACT comparison of the spec name against the tree name before the first
# ``|`` (never fuzzy). ``representative_name_overrides`` records a human-confirmed
# renamed sample. If a representative is genuinely absent from the tree, its clade
# must be listed in ``representative_absent_ok`` (so the cross-check is knowingly
# skipped for it); an unlisted absent representative is a hard error, not a guess.
# The config is version-controlled, so these entries are reviewed like code.
# ---------------------------------------------------------------------------


def validate_resolution_overrides(name_overrides, absent_ok):
    """Sanity-check the config-supplied override / absent-ok tables (fail loudly)."""
    if not isinstance(name_overrides, dict):
        raise ValueError("config 'representative_name_overrides' must be a mapping")
    for k, v in name_overrides.items():
        if not isinstance(v, str):
            raise ValueError(
                f"representative_name_overrides['{k}'] must be a string tree name")
    if not isinstance(absent_ok, (list, tuple)):
        raise ValueError(
            "config 'representative_absent_ok' must be a list of clade names")
    for c in absent_ok:
        if not isinstance(c, str):
            raise ValueError(
                "representative_absent_ok entries must be clade-name strings "
                f"(got {c!r})")


# ---------------------------------------------------------------------------
# Node-map intermediates exchanged between pipeline rules
#
# ``resolve_nodes.py`` writes these (one per anchoring mode); ``compare_and_build.py``
# reads them. The shared schema lives here so the two rules cannot drift apart.
# ---------------------------------------------------------------------------


def write_nodes(path, clade_to_node, via=None, unresolved=None, meta=None):
    """Write a small per-clade node-map YAML artifact.

    Schema::

        meta: {...}                 # optional free-form provenance
        nodes:
          <clade>: {node: <id>, via: <how it was resolved>}
        unresolved: [<clade>, ...]  # clades a mode could not anchor

    Kept human-readable and diffable so the two-mode node maps are easy to
    inspect between pipeline stages.
    """
    via = via or {}
    doc = {
        "meta": meta or {},
        "nodes": {c: {"node": clade_to_node[c], "via": via.get(c, "representative")}
                  for c in sorted(clade_to_node)},
        "unresolved": sorted(unresolved or []),
    }
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=True, default_flow_style=False)


def read_nodes(path):
    """Read a ``write_nodes`` artifact -> ``(clade_to_node, via, unresolved)``."""
    with open_maybe_gzip(path) as f:
        doc = yaml.safe_load(f) or {}
    nodes = doc.get("nodes", {}) or {}
    clade_to_node = {c: v["node"] for c, v in nodes.items()}
    via = {c: v.get("via", "representative") for c, v in nodes.items()}
    unresolved = list(doc.get("unresolved", []) or [])
    return clade_to_node, via, unresolved


# ---------------------------------------------------------------------------
# Config
#
# All run parameters come from a YAML config file (default: config.yaml). Paths
# are resolved relative to the working directory the command is run from (the
# repo root). This includes the human-approved representative_name_overrides /
# representative_absent_ok tables -- the config is version-controlled, so those
# entries are reviewed exactly like code.
# ---------------------------------------------------------------------------

# key -> default (keys with no default are required)
CONFIG_DEFAULTS = {
    "old_yaml": None,            # reference yaml for the old-scheme recovery check
    "clade_column": "Nextstrain_clade",
    "site_include_start": 266,   # masks outside [start, end] (the UTRs) are dropped
    "site_include_end": 29674,
    "matutils": "matUtils",
    "strict": False,             # exit nonzero if old masks are not fully recovered
    "log": None,
    # Second output: per-PANGO-lineage masks (companion to the Nextstrain output).
    "output_pango": "results/usher_masked_sites_pango.yaml",
    # Anchor cross-check: fail if a clade's founder anchor sits more than this many
    # mutations (patristic) from its representative anchor.
    "max_anchor_distance": 4,
    # Staged-pipeline intermediates / reports (defaults under results/; the
    # Snakemake rules write these and the downstream rules read them):
    "nodes_representative": "results/nodes_representative.yaml",  # Mode 1 node map
    "nodes_founder": "results/nodes_founder.yaml",               # Mode 2 node map
    "comparison_report": "results/anchoring_comparison.txt",     # anchor distance table
    "anchor_distances": "results/anchor_distances.yaml",         # machine-readable table
    "anchor_check_report": "results/anchor_check.txt",           # assertion result
    "recovery_report": "results/recovery_report.txt",            # old-scheme recovery
    "recovery_summary": "results/recovery_summary.txt",          # per-clade added-mask table
    # Human-approved resolution deviations (see validate_resolution_overrides):
    "representative_name_overrides": {},  # spec rep name -> exact tree name
    "representative_absent_ok": [],       # clades whose rep may be absent (skip check)
}
REQUIRED_CONFIG = ("branch_spec", "metadata", "output", "tree")


def load_config(path):
    """Load and validate the YAML config, returning a fully-defaulted dict."""
    with open_maybe_gzip(path) as f:
        raw = yaml.safe_load(f) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"Config {path} must be a YAML mapping")
    # Reject typos / stale keys up front rather than silently ignoring them.
    allowed = set(CONFIG_DEFAULTS) | set(REQUIRED_CONFIG)
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(
            f"Unknown config key(s) in {path}: {sorted(unknown)}. "
            f"Allowed keys: {sorted(allowed)}"
        )
    cfg = dict(CONFIG_DEFAULTS)
    cfg.update(raw)
    for key in REQUIRED_CONFIG:
        if not cfg.get(key):
            raise ValueError(f"Config {path} is missing required key '{key}'")
    return cfg
