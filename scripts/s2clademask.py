#!/usr/bin/env python3
"""Shared library for the s2-clade-mask generator.

This module holds every piece of logic used to turn the UCSC PANGO-lineage-keyed
``branchSpecificMask.yml`` into the Nextstrain-clade-keyed ``usher_masked_sites.yaml``
that ``mat-pcp`` consumes. It is imported by both:

* the staged Snakemake pipeline (``Snakefile`` + the thin ``scripts/*.py`` rule
  entrypoints ``extract_paths.py`` / ``resolve_nodes.py`` / ``compare_and_build.py``
  / ``check_recovery.py``), and
* the one-shot wrapper ``scripts/generate_usher_masked_sites.py``.

All the algorithm's write-up lives in ``README.md``; this module keeps the
functions and the design invariants (representative-primary resolution, exact
matching only, human-approved deviations, superset recovery, UTR drop) in one
place so the pipeline scripts stay thin and never re-implement logic.
"""

import os
import re
import subprocess
import gzip
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


# matUtils extract writes its output relative to its ``-d`` directory (default
# the CWD), prepending that directory even to an *absolute* ``--sample-paths`` /
# ``--clade-paths`` path -- so a bare absolute path lands in a bogus nested
# location (``$CWD/fh/.../out.tsv``) and the intended file is never written.
# Passing ``-d /`` together with an absolute out path makes the prepend collapse
# (``/`` + ``/abs`` -> ``//abs`` == ``/abs``), so the file lands exactly where we
# asked. Both extract wrappers must do this; forgetting it on one silently yields
# an empty output file.


def run_sample_paths(tree, out_path, matutils="matUtils", log_file=None):
    """Run ``matUtils extract --sample-paths`` on ``tree`` writing ``out_path``."""
    out_abs = os.path.abspath(out_path)
    cmd = [matutils, "extract", "-i", str(tree), "-d", "/", "--sample-paths", out_abs]
    log(f"Running: {' '.join(cmd)}", log_file)
    subprocess.run(cmd, check=True)
    return out_abs


def run_clade_paths(tree, out_path, matutils="matUtils", log_file=None):
    """Run ``matUtils extract --clade-paths`` on ``tree`` writing ``out_path``."""
    out_abs = os.path.abspath(out_path)
    cmd = [matutils, "extract", "-i", str(tree), "-d", "/", "--clade-paths", out_abs]
    log(f"Running: {' '.join(cmd)}", log_file)
    subprocess.run(cmd, check=True)
    return out_abs


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
# Human-approved resolution overrides (supplied via config)
#
# Representative samples are resolved to tree nodes by an EXACT comparison of the
# spec name against the tree sample name before the first ``|``. No fuzzy or
# substring matching is ever performed. When a spec representative does not match
# the tree verbatim, or is absent from the tree entirely, resolution requires a
# deliberate, human-reviewed entry in the config -- ``representative_name_overrides``
# (a renamed sample) or ``representative_fallbacks`` (a substituted node) -- so a
# maintainer signs off on every single deviation. The config is under version
# control, so these entries are reviewed exactly like code. A spec representative
# that is neither found nor listed is a hard error, not a silent guess.
# ---------------------------------------------------------------------------


def validate_resolution_overrides(name_overrides, fallbacks):
    """Sanity-check the config-supplied override/fallback tables (fail loudly)."""
    if not isinstance(name_overrides, dict):
        raise ValueError("config 'representative_name_overrides' must be a mapping")
    for k, v in name_overrides.items():
        if not isinstance(v, str):
            raise ValueError(
                f"representative_name_overrides['{k}'] must be a string tree name")
    if not isinstance(fallbacks, dict):
        raise ValueError("config 'representative_fallbacks' must be a mapping")
    for clade, fb in fallbacks.items():
        if not isinstance(fb, dict) or "via" not in fb:
            raise ValueError(
                f"representative_fallbacks['{clade}'] must be a mapping with a 'via' key")
        if fb["via"] not in ("founder", "sample"):
            raise ValueError(
                f"representative_fallbacks['{clade}'].via must be 'founder' or "
                f"'sample' (got {fb['via']!r})")
        if fb["via"] == "sample" and not fb.get("name"):
            raise ValueError(
                f"representative_fallbacks['{clade}'] with via=sample needs a 'name'")


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
            # Drop empty tokens: matUtils (usher >=0.6.6) emits a trailing space
            # after the last node, so a bare split yields a trailing '' that would
            # otherwise be picked up as the (empty) representative node.
            nodes = [tok for tok in path.split(" ") if tok]
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


def resolve_founder_nodes(spec, clade_founders):
    """Map every PANGO clade to its annotated founder (MRCA) node (Mode 2, pure).

    Returns ``(founder_node, no_founder)`` where ``founder_node`` maps each clade
    that has an annotated founder in the tree to that node id, and ``no_founder``
    is the sorted list of clades with no annotated founder (non-monophyletic, or
    simply not annotated in this tree) -- those cannot be founder-anchored.
    """
    founder_node = {}
    no_founder = []
    for clade in spec:
        node = clade_founders.get(clade)
        if node is None:
            no_founder.append(clade)
        else:
            founder_node[clade] = node
    return founder_node, sorted(no_founder)


def resolve_representative_only(spec, sample_paths_file, name_overrides,
                                log_file=None):
    """Mode 1 (string match): anchor each clade at its representative sample node.

    Steps 1-2 of the resolution policy, with **no** founder/sample fallback:

    1. **Representative sample (primary).** Resolve ``spec[clade]["representative"]``
       from the sample-paths file, matching the name *exactly* before the first
       ``|`` (upstream ``branchSpecificMask.py`` logic, respecting
       ``representativeBacktrack``). This anchors the mask at the same node
       ``matUtils mask`` would.
    2. **Approved name override.** If the spec name does not match the tree
       verbatim but a human has recorded the correct tree name in
       ``representative_name_overrides``, resolve using that name.

    Returns ``(clade_to_node, via, unresolved)``: the resolved nodes, how each was
    resolved (``"representative"`` / ``"name override"``), and the sorted list of
    clades whose representative is absent from the tree (candidates for an approved
    fallback in ``apply_fallbacks``). Run against the full **gisaid-and-public**
    MAT; the public-only tree is missing several GISAID-only representatives.
    """
    # Build the search name for each clade (approved override, or the spec name).
    search = {}
    overridden = {}
    for clade, entry in spec.items():
        rep = entry["representative"]
        name = name_overrides.get(rep, rep)
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
    unresolved = sorted(c for c in spec if c not in clade_to_node)

    for clade in sorted(overridden):
        spec_name, tree_name = overridden[clade]
        log(f"  name override (approved): {clade} '{spec_name}' -> '{tree_name}'",
            log_file)
    return clade_to_node, via, unresolved


def apply_fallbacks(spec, clade_to_node, via, clade_founders, sample_paths_file,
                    fallbacks, log_file=None):
    """Fill in clades unresolved by Mode 1 using human-approved config fallbacks.

    Steps 3-4 of the resolution policy. ``clade_to_node`` / ``via`` are Mode 1's
    output (mutated in place and also returned). For every clade in ``spec`` not
    already resolved:

    3. **Approved fallback** (``representative_fallbacks[clade]``): ``via: founder``
       uses the clade's annotated founder (MRCA) node -- the Mode-2 anchor -- and
       ``via: sample`` resolves an approved alternate sample from
       ``sample_paths_file``.
    4. Otherwise **fail loudly** naming the clade -- a maintainer must add an
       override or fallback; the script never guesses.
    """
    still_missing = []
    for clade in (c for c in spec if c not in clade_to_node):
        fb = fallbacks.get(clade)
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
                "config representative_fallbacks."
            )

    if still_missing:
        detail = "\n".join(
            f"    {c}: representative '{spec[c]['representative']}' not found in the tree"
            for c in sorted(still_missing)
        )
        raise ValueError(
            "Could not resolve a tree node for the following clade(s): their "
            "representative sample is absent from the tree and no approved "
            "override/fallback exists. Add an entry to the config's "
            "representative_name_overrides (if the sample is present under a "
            "different name) or representative_fallbacks (to approve the founder "
            f"node or an alternate sample):\n{detail}"
        )

    for clade in sorted(clade_to_node):
        log(f"  {clade} -> {clade_to_node[clade]} (via {via[clade]})", log_file)
    return clade_to_node, via


def resolve_clade_nodes(spec, clade_founders, sample_paths_file,
                        name_overrides, fallbacks, log_file=None):
    """Authoritative PANGO-clade -> tree-node resolver (Mode 1 + approved fallback).

    Composes ``resolve_representative_only`` (representative-primary, name
    overrides) with ``apply_fallbacks`` (approved per-clade founder/sample
    fallbacks; hard error otherwise). This is the single source of truth for the
    final anchoring, used by the one-shot ``generate_usher_masked_sites.py``; the
    staged pipeline reaches the same result by running the two steps as separate
    rules. Returns ``(clade_to_node, via)``.
    """
    clade_to_node, via, _ = resolve_representative_only(
        spec, sample_paths_file, name_overrides, log_file)
    return apply_fallbacks(spec, clade_to_node, via, clade_founders,
                           sample_paths_file, fallbacks, log_file)


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


def build_masks_for_nodes(spec, node_map, ns_full_values, clade_founders,
                          sample_paths_file, start, end, log_file=None):
    """Run the full over-mask assignment for one PANGO-clade->node anchoring.

    Convenience wrapper around ``assign_masks_by_founder`` -> ``build_clade_masks``
    -> ``restrict_to_include_range``. ``node_map`` maps ``pango_clade -> anchor
    node id`` (representatives, founders, or a mix). Returns the in-range
    ``{bare_ns_clade: {"sites", "mutations", "sources"}}`` map the generator emits.
    """
    ns_to_pango = assign_masks_by_founder(
        spec, node_map, ns_full_values, clade_founders, sample_paths_file, log_file)
    clade_masks = build_clade_masks(spec, ns_to_pango, log_file)
    restrict_to_include_range(clade_masks, start, end, log_file)
    return clade_masks


def diff_masks(rep_masks, founder_masks):
    """Diff two ``{ns_clade: {"sites","mutations"}}`` maps.

    Returns ``{ns_clade: {"sites_only_rep", "sites_only_founder", "muts_only_rep",
    "muts_only_founder"}}`` for the Nextstrain clades whose masks differ between
    the two anchorings.
    """
    diffs = {}
    for clade in sorted(set(rep_masks) | set(founder_masks)):
        r = rep_masks.get(clade, {})
        f = founder_masks.get(clade, {})
        r_sites, f_sites = r.get("sites", set()), f.get("sites", set())
        r_muts, f_muts = r.get("mutations", set()), f.get("mutations", set())
        entry = {
            "sites_only_rep": r_sites - f_sites,
            "sites_only_founder": f_sites - r_sites,
            "muts_only_rep": r_muts - f_muts,
            "muts_only_founder": f_muts - r_muts,
        }
        if any(entry.values()):
            diffs[clade] = entry
    return diffs


# ---------------------------------------------------------------------------
# Node-map intermediates exchanged between pipeline rules
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


def load_generated_masks(yaml_path):
    """Read a written ``usher_masked_sites.yaml`` into the in-memory mask map.

    Returns ``{bare_clade: {"sites": set[int], "mutations": set[str]}}`` -- the
    same shape ``build_clade_masks`` produces (minus ``sources``) -- so the
    standalone recovery rule can re-check the on-disk output without re-running
    the assignment. Blocks are keyed by their ``clades`` list, matching how
    ``load_old_sites`` reads the old scheme.
    """
    with open_maybe_gzip(yaml_path) as f:
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
# Config
#
# All run parameters come from a YAML config file (default: config.yaml). Paths
# are resolved relative to the working directory the command is run from (the
# repo root). This includes the human-approved representative_name_overrides /
# representative_fallbacks tables -- the config is version-controlled, so those
# entries are reviewed exactly like code.
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
    # Staged-pipeline intermediates / reports (defaults under workdir; the
    # Snakemake rules write these, and compare_and_build/check_recovery read them):
    "nodes_representative": "results/nodes_representative.yaml",  # Mode 1 node map
    "nodes_founder": "results/nodes_founder.yaml",               # Mode 2 node map
    "comparison_report": "results/anchoring_comparison.txt",     # rep-vs-founder diff
    "recovery_report": "results/recovery_report.txt",            # old-scheme recovery
    # Human-approved resolution deviations (see validate_resolution_overrides):
    "representative_name_overrides": {},  # spec rep name -> exact tree name
    "representative_fallbacks": {},       # clade -> {via: founder|sample, ...}
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
