#!/usr/bin/env python3
"""Generate the Nextstrain-clade-keyed ``usher_masked_sites.yaml`` from the new
UCSC PANGO-lineage-keyed ``branchSpecificMask.yml``.

Background
----------
The UShER tree-building pipeline replaced its old, hand-maintained per-clade
masking (``maskDelta.sh``) with a data-driven, tree-topology-based spec
(``branchSpecificMask.yml`` + ``branchSpecificMask.py``); see
https://github.com/yatisht/usher/issues/324 . In the new scheme each mask is
keyed by a **PANGO lineage**, anchored to a **representative sample**, and
applied by ``matUtils mask`` to that sample's node **and its whole subtree** --
so inheritance down the tree is implicit via topology.

``mat-pcp`` works at the level of **Nextstrain clades** and masks by exact clade
name (``data/usher_masked_sites.yaml``, consumed in ``count_mutations.py`` and
``aggregate_mutations_to_exclude.py``). This script bridges the two: it resolves
each PANGO clade's representative to a tree node (reusing the exact logic of the
reference ``branchSpecificMask.py``), uses the tree to determine which Nextstrain
clades fall under each representative node, and emits a flat, fully-accumulated
``Nextstrain_clade -> {sites, mutations}`` map.

Approach (tree topology + over-mask)
------------------------------------
* Each PANGO clade is anchored at its **representative sample's node** -- the same
  node ``matUtils mask`` applies the mask to upstream. The representative is
  resolved from ``matUtils extract --sample-paths`` by an *exact* match on the
  sample name before ``|`` (respecting ``representativeBacktrack``), never a fuzzy
  or substring match. When a representative does not match the tree verbatim, or
  is absent from it, resolution requires a deliberate, hard-coded, human-reviewed
  entry (``REPRESENTATIVE_NAME_OVERRIDES`` / ``REPRESENTATIVE_FALLBACKS``); an
  unlisted, unresolved representative is a hard error. Run against a tree that
  contains the representatives -- the full **gisaid-and-public** MAT, since several
  representatives are GISAID-only and absent from the public-only tree.
* Nextstrain clades are anchored at their annotated **founder (MRCA) nodes** from
  ``matUtils extract --clade-paths``. The set of Nextstrain clades to emit comes
  from the metadata ``clade_column``; each full value (e.g. ``"22A (BA.4)"``) keys
  into ``--clade-paths`` and is reduced to its bare token (``"22A"``) for output,
  exactly as ``get_clades`` does.
* Over-mask rule: a Nextstrain clade ``C`` inherits PANGO clade ``P``'s masks when
  ``P``'s representative node and ``C``'s founder lie on one root->tip lineage --
  either ``P`` is an ancestor of ``C`` (all of ``C`` descends from ``P``; exact) or
  ``P`` is nested inside ``C`` (finer PANGO clade, e.g. ``BQ.1.1`` inside ``22E``;
  widened to the whole Nextstrain clade since mat-pcp cannot mask below Nextstrain
  granularity). A consequence is that trunk clades like ``21L`` (BA.2) accumulate
  their descendant sub-clades' masks. Non-monophyletic clades (``recombinant``)
  have no founder node and are left unmasked. Contributions are logged per clade.
* ``matUtils --sample-paths`` is also used to read each anchor node's ancestor
  chain from the tree.
* ``ranges`` and ``sites`` become masked **sites**. ``reversions`` become masked
  **mutations** (directional ``refPOSalt`` strings, matching mat-pcp's
  ``nt_mutation``); they are *not* expanded to whole-site masks. ``exclusions``
  (per-sample un-masking of a site) cannot be expressed at Nextstrain
  granularity and are logged and skipped.

Recovering the old scheme
-------------------------
The generated file must be a **superset** of the current hand-made one: for
every ``(clade, site)`` masked before, the new file must still mask it. The
script checks ``old[clade] ⊆ new[clade]`` and reports any missing pairs. In
the ``strict`` config option makes a shortfall a hard error to investigate;
``union_old`` force-includes the old masks (guaranteeing recovery) while logging
the pairs that came only from the old file.

Resources
---------
This is an offline generator meant to run occasionally against the **full
gisaid-and-public** MAT (all clades must be present/annotated -- do not run it
against a filtered per-run tree; several representative samples are GISAID-only
and absent from the public tree). ``matUtils extract --sample-paths`` on an
~8M-sample tree produces a large paths file (tens of GB); the generator only reads
it until every anchor node's ancestor chain has been seen (early exit), so memory
stays modest. Point the ``workdir`` config option at a disk with room, or pass
precomputed ``sample_paths`` / ``clade_paths`` to skip matUtils.

Usage
-----
All run parameters come from a YAML config (default ``config.yaml``); run from the
repo root so relative paths resolve:

    python scripts/generate_usher_masked_sites.py --config config.yaml
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

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
# 1. Parse the new UCSC spec
# ---------------------------------------------------------------------------

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

    version = raw.pop("version", None)
    if version is None:
        raise ValueError(f"No 'version' key found in {path}")

    spec = {}
    for clade, entry in raw.items():
        rep = entry.get("representative")
        if not rep:
            raise ValueError(f"No 'representative' for clade {clade} in {path}")

        sites = set()
        for lo, hi in entry.get("ranges", []) or []:
            sites.update(range(int(lo), int(hi) + 1))
        for s in entry.get("sites", []) or []:
            sites.add(int(s))

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
# 2. Sample paths (root -> tip node IDs) from the MAT
# ---------------------------------------------------------------------------


def run_sample_paths(tree, out_path, matutils="matUtils", log_file=None):
    """Run ``matUtils extract --sample-paths`` on ``tree`` writing ``out_path``."""
    cmd = [matutils, "extract", "-i", str(tree), "-d", "/", "--sample-paths", str(out_path)]
    log(f"Running: {' '.join(cmd)}", log_file)
    subprocess.run(cmd, check=True)
    return out_path


def run_clade_paths(tree, out_path, matutils="matUtils", log_file=None):
    """Run ``matUtils extract --clade-paths`` on ``tree`` writing ``out_path``."""
    cmd = [matutils, "extract", "-i", str(tree), "--clade-paths", str(out_path)]
    log(f"Running: {' '.join(cmd)}", log_file)
    subprocess.run(cmd, check=True)
    return out_path


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
            if len(parts) >= 2 and parts[0] and parts[1]:
                founders.setdefault(parts[0], parts[1])
    return founders


# ---------------------------------------------------------------------------
# Human-approved resolution overrides
#
# Representative samples are resolved to tree nodes by an EXACT comparison of the
# spec name against the tree sample name before the first ``|``. No fuzzy or
# substring matching is ever performed. When a spec representative does not match
# the tree verbatim, or is absent from the tree entirely, resolution requires a
# deliberate, hard-coded entry below -- so a maintainer signs off on every single
# deviation (a renamed sample or a substituted node). A spec representative that
# is neither found nor listed here is a hard error, not a silent guess.
# ---------------------------------------------------------------------------

# Spec representative name -> exact tree sample name (the field before ``|``).
# Add an entry only after confirming by hand that the two names denote the same
# sample. Verified against the 2026-04-22 gisaid-and-public MAT, where these
# differ from the spec only by country abbreviation.
REPRESENTATIVE_NAME_OVERRIDES = {
    "IND/GBRC714b/2021": "India/GJ-GBRC714b/2021",                        # clade B.1.617.2
    "BGD/icddrb_TND_06_1053/2022": "Bangladesh/icddrb_TND_06_1053/2022",  # clade XBB
}

# PANGO clade -> approved fallback used only when the clade's representative
# sample is absent from the tree. ``{"via": "founder"}`` uses the clade's
# annotated founder (MRCA) node; ``{"via": "sample", "name": ..., "backtrack": n}``
# resolves an approved alternate sample. Each entry is an explicit human decision
# to resolve that specific clade some other way.
REPRESENTATIVE_FALLBACKS = {
    # FRA/IHUCOVID-005193-N1/2021 is not present in the 2026-04-22
    # gisaid-and-public MAT. P.1's only mask is the deletion [11288, 11296], which
    # is also masked via B.1.1.7 and B.1.351; the annotated 20J (Gamma) founder
    # anchors it faithfully.
    "P.1": {"via": "founder"},
}


def resolve_representative_nodes(clade_search, sample_paths_file, log_file=None):
    """Resolve each clade's representative sample to a tree node ID.

    ``clade_search`` maps ``clade -> {"name": tree_name, "backtrack": int}``, where
    ``tree_name`` is the exact sample name to match **before the first** ``|``.
    Mirrors ``branchSpecificMask.py``'s ``getRepresentativeNodes``: take the last
    node in the sample's path, back up one if it is the sample's own leaf, then
    back up ``backtrack`` more. Matching is strict equality on the pre-``|`` field
    -- never a substring match. Returns ``{clade: node_id}`` for the samples that
    were actually found in the tree (clades whose sample is absent are omitted).
    """
    # tree sample name -> list of (clade, backtrack); a name may anchor >1 clade.
    wanted = defaultdict(list)
    for clade, info in clade_search.items():
        wanted[info["name"]].append((clade, info["backtrack"]))

    clade_to_node = {}
    remaining = set(wanted)
    with open_maybe_gzip(sample_paths_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if "\t" not in line:
                continue
            full_name, path = line.split("\t", 1)
            name = full_name.split("|")[0]
            if name not in remaining:
                continue
            nodes = path.split(" ")
            for clade, backtrack in wanted[name]:
                idx = -1
                if nodes[idx].startswith(name):
                    idx -= 1
                idx -= backtrack
                if -idx > len(nodes):
                    raise ValueError(
                        f"Backtrack for representative '{name}' (clade {clade}) "
                        f"walks off the root (path has {len(nodes)} nodes)"
                    )
                clade_to_node[clade] = nodes[idx].split(":", 1)[0]
            remaining.discard(name)
            if not remaining:
                break

    return clade_to_node


def resolve_clade_nodes(spec, clade_founders, sample_paths_file, log_file=None):
    """Map each PANGO clade to a tree node ID via its representative sample.

    Resolution is deliberately explicit and auditable; nothing is matched fuzzily
    and no substitution happens without a hard-coded, human-reviewed entry:

    1. **Representative sample (primary).** Resolve ``spec[clade]["representative"]``
       from the sample-paths file, matching the name *exactly* before the first
       ``|`` (upstream ``branchSpecificMask.py`` logic, respecting
       ``representativeBacktrack``). This anchors the mask at the same node
       ``matUtils mask`` would.
    2. **Approved name override.** If the spec name does not match the tree
       verbatim but a human has recorded the correct tree name in
       ``REPRESENTATIVE_NAME_OVERRIDES``, resolve using that name.
    3. **Approved fallback.** If the representative is genuinely absent from the
       tree and a human has approved an alternative in ``REPRESENTATIVE_FALLBACKS``,
       apply it (the clade's annotated founder node, or an approved alternate
       sample).
    4. Otherwise **fail loudly** -- a maintainer must add an override or fallback.

    Run against a tree that actually contains the representative samples (the full
    **gisaid-and-public** MAT); the public-only tree is missing several
    GISAID-only representatives.
    """
    # Build the search name for each clade (approved override, or the spec name).
    search = {}
    overridden = {}
    for clade, entry in spec.items():
        rep = entry["representative"]
        name = REPRESENTATIVE_NAME_OVERRIDES.get(rep, rep)
        if name != rep:
            overridden[clade] = (rep, name)
        search[clade] = {"name": name, "backtrack": entry["backtrack"]}

    resolved = resolve_representative_nodes(search, sample_paths_file, log_file)

    clade_to_node = {}
    via = {}
    for clade in spec:
        if clade in resolved:
            clade_to_node[clade] = resolved[clade]
            via[clade] = "name override" if clade in overridden else "representative"

    # Representatives that were not found require an approved fallback.
    still_missing = []
    for clade in (c for c in spec if c not in clade_to_node):
        fb = REPRESENTATIVE_FALLBACKS.get(clade)
        if fb is None:
            still_missing.append(clade)
            continue
        kind = fb.get("via")
        if kind == "founder":
            node = clade_founders.get(clade)
            if node is None:
                raise ValueError(
                    f"Approved 'founder' fallback for clade {clade}, but the clade "
                    "has no annotated founder node in the tree."
                )
            clade_to_node[clade] = node
            via[clade] = "approved founder fallback"
        elif kind == "sample":
            alt = fb["name"]
            alt_res = resolve_representative_nodes(
                {clade: {"name": alt, "backtrack": fb.get("backtrack", 0)}},
                sample_paths_file, log_file)
            if clade not in alt_res:
                raise ValueError(
                    f"Approved alternate sample '{alt}' for clade {clade} was not "
                    "found in the tree either."
                )
            clade_to_node[clade] = alt_res[clade]
            via[clade] = f"approved alternate sample ({alt})"
        else:
            raise ValueError(
                f"Unknown fallback kind {kind!r} for clade {clade} in "
                "REPRESENTATIVE_FALLBACKS."
            )

    if still_missing:
        detail = "\n".join(
            f"    {c}: representative '{spec[c]['representative']}' not found in the tree"
            for c in sorted(still_missing)
        )
        raise ValueError(
            "Could not resolve a tree node for the following clade(s): their "
            "representative sample is absent from the tree and no approved "
            "override/fallback exists. Add an entry to REPRESENTATIVE_NAME_OVERRIDES "
            "(if the sample is present under a different name) or "
            "REPRESENTATIVE_FALLBACKS (to approve the founder node or an alternate "
            f"sample):\n{detail}"
        )

    for clade in sorted(overridden):
        spec_name, tree_name = overridden[clade]
        log(f"  name override (approved): {clade} '{spec_name}' -> '{tree_name}'",
            log_file)
    for clade in sorted(clade_to_node):
        log(f"  {clade} -> {clade_to_node[clade]} (via {via[clade]})", log_file)
    return clade_to_node


# ---------------------------------------------------------------------------
# 3. Metadata: strain -> bare Nextstrain clade
# ---------------------------------------------------------------------------


def bare_clade(value):
    """Reduce a clade field to its bare token, as ``get_clades`` does.

    ``"22A (BA.4)"`` -> ``"22A"``; ``"recombinant"`` -> ``"recombinant"``.
    """
    return value.split(" ", 1)[0]


def load_nextstrain_clade_values(metadata, clade_column, log_file=None):
    """Return the set of distinct full ``clade_column`` values in the metadata.

    Values are the full annotations (e.g. ``"22A (BA.4)"``); the bare token
    (``"22A"``) is used as the output key, matching the ``get_clades`` rule.
    """
    values = set()
    with open_maybe_gzip(metadata) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            clade_i = header.index(clade_column)
        except ValueError:
            raise ValueError(f"Metadata is missing the '{clade_column}' column")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > clade_i and fields[clade_i]:
                values.add(fields[clade_i])
    log(f"Found {len(values)} distinct '{clade_column}' values in metadata",
        log_file)
    return values


# ---------------------------------------------------------------------------
# 4. Assign PANGO masks to Nextstrain clades via tree topology (over-mask)
# ---------------------------------------------------------------------------


def capture_ancestor_chains(sample_paths_file, target_nodes, log_file=None):
    """Return ``{node: [ancestor node ids from root down to and incl. node]}``.

    Ancestry is read from the tree topology: for each target founder node, the
    first sample-path that passes through it yields the ordered chain of nodes
    from the root down to (and including) that node. Reads only until every
    target is found (early exit), so this is cheap even on a huge paths file.
    """
    remaining = set(target_nodes)
    chains = {}
    with open_maybe_gzip(sample_paths_file) as f:
        for line in f:
            tab = line.find("\t")
            if tab < 0:
                continue
            ids = [tok.split(":", 1)[0] for tok in line[tab + 1:].split(" ") if tok]
            for node in remaining.intersection(ids):
                chains[node] = ids[: ids.index(node) + 1]
                remaining.discard(node)
            if not remaining:
                break
    if remaining:
        log(f"  WARNING: {len(remaining)} founder node(s) never appeared in any "
            f"sample path: {sorted(remaining)}", log_file)
    return chains


def assign_masks_by_founder(spec, pango_node, ns_full_values, clade_founders,
                            sample_paths_file, log_file=None):
    """Return ``{bare_nextstrain_clade: set(pango_clade_names)}``.

    Each PANGO clade ``P`` is anchored at its **representative node** ``r(P)`` (in
    ``pango_node``) -- the same node ``matUtils mask`` applies the mask to upstream
    -- while each Nextstrain clade ``C`` is anchored at its annotated **founder**
    (MRCA) node. ``C`` inherits ``P``'s masks when either:

    (a) ``P`` is an ancestor of ``C`` (``r(P)`` lies on the root->founder(C)
        path) -- every sample of ``C`` descends from ``r(P)``, so the mask applies
        exactly; or
    (b) ``C`` is the *deepest* Nextstrain clade enclosing ``P`` (``C``'s founder is
        the last Nextstrain founder on the root->``r(P)`` path) -- ``P``'s "stem"
        samples, i.e. those not relabeled into a finer Nextstrain clade, carry the
        Nextstrain label ``C``, so ``P``'s mask is widened to all of ``C``
        (over-mask coarsening).

    Working from the representative node (not raw member samples) keeps this free
    of the label noise -- recombinants, misannotations -- that a per-sample
    membership rule would let bleed masks across unrelated clades. Rule (b) also
    ensures a masked PANGO clade is attributed to the clade its samples actually
    carry (e.g. ``KP.3.1.1`` -> the ``24E`` clade, not its parent ``24A``), rather
    than to every ancestor.
    """
    ns_founder = {}
    skipped = []
    for full in ns_full_values:
        node = clade_founders.get(full)
        if node is None:
            skipped.append(full)
        else:
            ns_founder[full] = node
    if skipped:
        log(f"  {len(skipped)} Nextstrain clade(s) have no tree founder node "
            f"(not monophyletic, e.g. recombinant) -- left unmasked: "
            f"{sorted(skipped)}", log_file)

    # node -> Nextstrain clade full-name(s) with that founder (usually one).
    nsnode_to_full = defaultdict(list)
    for full, node in ns_founder.items():
        nsnode_to_full[node].append(full)
    ns_nodes = set(ns_founder.values())

    targets = set(pango_node.values()) | ns_nodes
    chains = capture_ancestor_chains(sample_paths_file, targets, log_file)
    chain_sets = {n: set(c) for n, c in chains.items()}

    ns_to_pango = defaultdict(set)

    # (a) P ancestor-or-equal of C: founder(P) on the root->founder(C) chain.
    for full, fc in ns_founder.items():
        bare = bare_clade(full)
        anc_c = chain_sets.get(fc, set())
        for p, fp in pango_node.items():
            if fp in anc_c:
                ns_to_pango[bare].add(p)

    # (b) attribute P to the deepest Nextstrain clade enclosing founder(P).
    for p, fp in pango_node.items():
        deepest = None
        for node in chains.get(fp, []):  # root -> founder(P) order
            if node in ns_nodes:
                deepest = node
        if deepest is None:
            log(f"  note: PANGO clade {p} has no enclosing Nextstrain clade in "
                "the tree; its stem samples are unlabeled and left unmasked",
                log_file)
            continue
        for full in nsnode_to_full[deepest]:
            ns_to_pango[bare_clade(full)].add(p)

    return ns_to_pango


def build_clade_masks(spec, ns_to_pango, log_file=None):
    """Accumulate per-Nextstrain-clade sites and mutations from assigned PANGO clades."""
    clade_masks = {}
    for ns_clade, pango_clades in sorted(ns_to_pango.items()):
        sites = set()
        mutations = set()
        for p in pango_clades:
            sites.update(spec[p]["sites"])
            mutations.update(spec[p]["mutations"])
        clade_masks[ns_clade] = {
            "sites": sites,
            "mutations": mutations,
            "sources": sorted(pango_clades),
        }
        log(f"  {ns_clade}: {len(sites)} sites, {len(mutations)} mutations "
            f"from {sorted(pango_clades)}", log_file)
    return clade_masks


def restrict_to_include_range(clade_masks, start, end, log_file=None):
    """Drop sites and reversion mutations outside ``[start, end]`` in place.

    The pipeline only ever considers positions within its ``site_include_range``
    (everything else -- notably the 5'/3' UTRs -- is discarded by the range
    filter at runtime). Masking those positions here is redundant, so for
    simplicity the generated file lists only in-range masks. Returns the number
    of dropped (clade, site) and (clade, mutation) entries.
    """
    dropped_sites = dropped_muts = 0
    for block in clade_masks.values():
        in_sites = {s for s in block["sites"] if start <= s <= end}
        dropped_sites += len(block["sites"]) - len(in_sites)
        block["sites"] = in_sites
        in_muts = {m for m in block["mutations"] if start <= int(m[1:-1]) <= end}
        dropped_muts += len(block["mutations"]) - len(in_muts)
        block["mutations"] = in_muts
    log(f"Restricted to site_include_range [{start}, {end}]: dropped "
        f"{dropped_sites} UTR/out-of-range site entries and {dropped_muts} "
        f"mutation entries.", log_file)
    return dropped_sites, dropped_muts


# ---------------------------------------------------------------------------
# 5. Old-scheme recovery check
# ---------------------------------------------------------------------------


def load_old_sites(old_yaml):
    """Return ``{clade: set(sites)}`` from an existing usher_masked_sites.yaml."""
    with open_maybe_gzip(old_yaml) as f:
        old = yaml.safe_load(f)
    result = defaultdict(set)
    for _, block in old.items():
        clades = block.get("clades", []) or []
        sites = set(block.get("sites", []) or [])
        for c in clades:
            result[c].update(sites)
    return result


def check_recovers_old(old_sites, clade_masks, union_old, include_range=None,
                       log_file=None):
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
        log(f"OLD-SCHEME RECOVERY: {n} old (clade, site) mask(s) fall outside the "
            f"site_include_range [{lo}, {hi}] and are intentionally NOT recovered "
            "(UTR sites dropped; discarded by the range filter at runtime).",
            log_file)

    if downgraded:
        n = sum(len(v) for v in downgraded.values())
        log(f"OLD-SCHEME RECOVERY: {n} old whole-site mask(s) are now masked at "
            "the mutation level (intentional site->mutation downgrade):", log_file)
        for clade in sorted(downgraded):
            log(f"  {clade}: {sorted(downgraded[clade])}", log_file)

    if truly_missing:
        log("OLD-SCHEME RECOVERY: the following (clade, site) masks are present "
            "in the old file but NOT recovered (neither site nor mutation):",
            log_file)
        for clade in sorted(truly_missing):
            log(f"  {clade}: {sorted(truly_missing[clade])}", log_file)
        if union_old:
            log("--union-old set: force-including the unrecovered old masks.",
                log_file)
            for clade, gap in truly_missing.items():
                clade_masks.setdefault(
                    clade, {"sites": set(), "mutations": set(), "sources": []}
                )
                clade_masks[clade]["sites"].update(gap)
                if "old:usher_masked_sites.yaml" not in clade_masks[clade]["sources"]:
                    clade_masks[clade]["sources"].append("old:usher_masked_sites.yaml")
    else:
        log("OLD-SCHEME RECOVERY: OK -- every old (clade, site) mask is recovered "
            "as a site or mutation mask.", log_file)
    return truly_missing


# ---------------------------------------------------------------------------
# 6. Emit usher_masked_sites.yaml
# ---------------------------------------------------------------------------


def write_output(out_path, clade_masks, spec, provenance, dropped_exclusions):
    """Write the Nextstrain-clade-keyed masking YAML, one block per clade."""
    lines = []
    lines.append("# ============================================================"
                 "=================")
    lines.append("# GENERATED FILE -- DO NOT EDIT BY HAND.")
    lines.append("#")
    lines.append("# Produced by the s2-clade-mask generator generate_usher_masked_sites.py")
    lines.append("# from the PANGO-lineage-keyed UCSC spec data/branchSpecificMask.yml, resolved")
    lines.append("# against a full gisaid-and-public MAT + metadata. Masks are accumulated down the")
    lines.append("# tree (topology-based inheritance) and coarsened to Nextstrain clades")
    lines.append("# (over-mask policy for PANGO clades finer than a Nextstrain clade).")
    lines.append("# Sites and reversion mutations outside the pipeline's")
    lines.append("# site_include_range (the UTRs) are omitted for simplicity -- they")
    lines.append("# are discarded by the range filter at runtime regardless.")
    lines.append("#")
    for line in provenance:
        lines.append(f"# {line}")
    if dropped_exclusions:
        lines.append("#")
        lines.append("# NOTE: the following per-sample 'exclusions' from the source spec")
        lines.append("# cannot be represented at Nextstrain-clade granularity and were")
        lines.append("# dropped (their sites remain masked clade-wide):")
        for clade, exs in dropped_exclusions.items():
            for ex in exs:
                lines.append(f"#   {clade}: site {ex['site']} "
                             f"(rep {ex['representative']})")
    lines.append("# ============================================================"
                 "=================")
    lines.append("")

    for i, clade in enumerate(sorted(clade_masks), start=1):
        block = clade_masks[clade]
        sites = sorted(block["sites"])
        mutations = sorted(block["mutations"])
        if not sites and not mutations:
            continue
        lines.append(f"mask_{clade}:")
        lines.append(f"  # accumulated from PANGO clades: {block['sources']}")
        lines.append("  clades:")
        lines.append(f"    - {clade}")
        if sites:
            lines.append(f"  sites: {sites}")
        if mutations:
            lines.append(f"  mutations: {mutations}")
        lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Config
#
# All run parameters come from a YAML config file (default: config.yaml). Paths
# are resolved relative to the working directory the command is run from (the
# repo root). The human-approved REPRESENTATIVE_NAME_OVERRIDES / _FALLBACKS
# tables are deliberately NOT config -- they are reviewed code decisions.
# ---------------------------------------------------------------------------

# key -> default (keys with no default are required)
CONFIG_DEFAULTS = {
    "old_yaml": None,            # reference yaml for the old-scheme recovery check
    "clade_column": "Nextstrain_clade",
    "site_include_start": 266,   # masks outside [start, end] (the UTRs) are dropped
    "site_include_end": 29674,
    "sample_paths": None,        # precomputed matUtils --sample-paths (skips matUtils)
    "clade_paths": None,         # precomputed matUtils --clade-paths (skips matUtils)
    "tree": None,                # full MAT .pb(.gz); required unless both *_paths given
    "matutils": "matUtils",
    "workdir": "results",        # dir for the (large) temporary paths files
    "union_old": False,          # force-include any unrecovered old masks
    "strict": False,             # exit nonzero if old masks are not fully recovered
    "log": None,
}
REQUIRED_CONFIG = ("branch_spec", "metadata", "output")


def load_config(path):
    """Load and validate the YAML config, returning a fully-defaulted dict."""
    with open_maybe_gzip(path) as f:
        raw = yaml.safe_load(f) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"Config {path} must be a YAML mapping")
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml",
                   help="Path to the YAML config (default: config.yaml). See the "
                        "repo's config.yaml for the schema.")
    cli = p.parse_args(argv)

    args = argparse.Namespace(**load_config(cli.config))
    log_file = args.log
    log(f"Loaded config: {cli.config}", log_file)

    # Ensure the scratch dir (for the large paths files) and output dir exist.
    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    log(f"Parsing branch spec: {args.branch_spec}", log_file)
    spec, version = parse_branch_spec(args.branch_spec)
    log(f"  {len(spec)} PANGO clades, spec version {version}", log_file)

    # Obtain sample paths (root->tip node IDs for every sample; used for the
    # membership scan and as a fallback representative-sample resolver).
    tmp_paths = None
    sample_paths = args.sample_paths
    if sample_paths is None:
        if not args.tree:
            raise ValueError("config needs 'tree' (or a precomputed 'sample_paths')")
        fd, tmp_paths = tempfile.mkstemp(prefix="sample_paths_", suffix=".tsv",
                                         dir=args.workdir)
        os.close(fd)
        sample_paths = run_sample_paths(args.tree, tmp_paths, args.matutils, log_file)

    # Obtain clade paths (annotated founder node per clade; primary node source).
    tmp_clade_paths = None
    clade_paths = args.clade_paths
    if clade_paths is None:
        if not args.tree:
            raise ValueError("config needs 'tree' (or a precomputed 'clade_paths')")
        fd, tmp_clade_paths = tempfile.mkstemp(prefix="clade_paths_", suffix=".tsv",
                                               dir=args.workdir)
        os.close(fd)
        clade_paths = run_clade_paths(args.tree, tmp_clade_paths, args.matutils,
                                      log_file)
    clade_founders = parse_clade_paths(clade_paths)
    log(f"Parsed {len(clade_founders):,} clade annotations from {clade_paths}",
        log_file)

    try:
        log("Resolving PANGO clade nodes (representative sample primary; "
            "approved overrides/fallbacks only)...", log_file)
        clade_to_node = resolve_clade_nodes(spec, clade_founders, sample_paths,
                                            log_file)

        log(f"Loading Nextstrain clades from {args.metadata}...", log_file)
        ns_full_values = load_nextstrain_clade_values(args.metadata,
                                                      args.clade_column, log_file)

        log("Assigning PANGO masks to Nextstrain clades "
            "(ancestor-or-nested over-mask; PANGO clades anchored at their "
            "representative node)...", log_file)
        ns_to_pango = assign_masks_by_founder(spec, clade_to_node, ns_full_values,
                                              clade_founders, sample_paths, log_file)
        clade_masks = build_clade_masks(spec, ns_to_pango, log_file)
        restrict_to_include_range(clade_masks, args.site_include_start,
                                  args.site_include_end, log_file)
    finally:
        if tmp_paths and os.path.exists(tmp_paths):
            os.remove(tmp_paths)
        if tmp_clade_paths and os.path.exists(tmp_clade_paths):
            os.remove(tmp_clade_paths)

    # Report which PANGO clades contribute masks to each Nextstrain clade. When
    # a Nextstrain clade lists PANGO clades finer than itself, its mask has been
    # widened to the whole clade (over-mask coarsening).
    log("Mask sources per Nextstrain clade (contributing PANGO clades):", log_file)
    for ns_clade in sorted(clade_masks):
        sources = [s for s in clade_masks[ns_clade]["sources"]
                   if not s.startswith("old:")]
        if sources:
            log(f"  {ns_clade} <- {sources}", log_file)

    dropped_exclusions = {c: e["exclusions"] for c, e in spec.items() if e["exclusions"]}

    missing = {}
    if args.old_yaml:
        old_sites = load_old_sites(args.old_yaml)
        missing = check_recovers_old(
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
    write_output(args.output, clade_masks, spec, provenance, dropped_exclusions)
    log(f"Wrote {args.output}", log_file)

    if missing and args.strict and not args.union_old:
        log("ERROR: --strict and old masks not fully recovered (see above).",
            log_file)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
