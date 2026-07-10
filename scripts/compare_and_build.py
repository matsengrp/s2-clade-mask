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

The over-mask assignment, anchor-distance, and YAML-emission logic all live in this
script (they are used by no other rule); the shared parsers / config / node-map
schema come from ``s2clademask``.

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
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


# ---------------------------------------------------------------------------
# Metadata: full clade value -> bare Nextstrain clade
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
    with lib.open_maybe_gzip(metadata) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            clade_i = header.index(clade_column)
        except ValueError:
            raise ValueError(f"Metadata is missing the '{clade_column}' column")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > clade_i and fields[clade_i]:
                values.add(fields[clade_i])
    lib.log(f"Found {len(values)} distinct '{clade_column}' values in metadata",
            log_file)
    return values


# ---------------------------------------------------------------------------
# Ancestor chains from the sample-paths file (tree topology)
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
    with lib.open_maybe_gzip(sample_paths_file) as f:
        for line in f:
            tab = line.find("\t")
            if tab < 0:
                continue
            ids = [tok.split(":", 1)[0] for tok in line[tab + 1:].split(" ") if tok]
            # A single path may pass through several targets at once; record each.
            for node in remaining.intersection(ids):
                chains[node] = ids[: ids.index(node) + 1]
                remaining.discard(node)
            if not remaining:
                break
    if remaining:
        lib.log(f"  WARNING: {len(remaining)} founder node(s) never appeared in any "
                f"sample path: {sorted(remaining)}", log_file)
    return chains


def capture_token_chains(sample_paths_file, target_nodes, log_file=None):
    """Like ``capture_ancestor_chains`` but keeps each node's mutation tokens.

    Returns ``{node: [tokens]}`` where each token is ``node_id:MUT,MUT`` (or bare
    ``node_id``) along the root->node path, so branch mutation counts are available.
    Strip the ``:MUT`` payload from each token to recover the id-only chains
    ``capture_ancestor_chains`` returns, so one token scan can feed both the anchor
    cross-check and the over-mask assignment.
    """
    remaining = set(target_nodes)
    chains = {}
    with lib.open_maybe_gzip(sample_paths_file) as f:
        for line in f:
            tab = line.find("\t")
            if tab < 0:
                continue
            toks = [t for t in line[tab + 1:].split(" ") if t]
            ids = [t.split(":", 1)[0] for t in toks]
            for node in remaining.intersection(ids):
                chains[node] = toks[: ids.index(node) + 1]
                remaining.discard(node)
            if not remaining:
                break
    if remaining:
        lib.log(f"  WARNING: {len(remaining)} anchor node(s) never appeared in any "
                f"sample path: {sorted(remaining)}", log_file)
    return chains


def _token_muts(token):
    """Mutation strings carried by one ``node_id:MUT,MUT`` sample-path token."""
    if ":" not in token:
        return []
    payload = token.split(":", 1)[1]
    return [m for m in payload.split(",") if m]


# ---------------------------------------------------------------------------
# Over-mask assignment: attribute PANGO masks to target clades via topology
# ---------------------------------------------------------------------------


def assign_masks_by_founder(spec, pango_node, ns_full_values, clade_founders,
                            sample_paths_file, log_file=None, chains=None):
    """Return ``{bare_target_clade: set(pango_clade_names)}``.

    ``pango_node`` maps each source PANGO clade to its anchor node; ``ns_full_values``
    are the *target* clade labels to attribute masks to (Nextstrain full values like
    ``"22A (BA.4)"`` for the Nextstrain output, or the PANGO clade names themselves
    for the PANGO-level output). Each target's founder (MRCA) node comes from
    ``clade_founders``. Pass a precomputed ``chains`` dict (from
    ``capture_ancestor_chains``, covering every anchor + target founder node) to
    reuse a single sample-paths scan across both outputs.

    Each PANGO clade ``P`` is anchored at node ``r(P)`` (in ``pango_node``); each
    target clade ``C`` is anchored at its founder node. ``C`` inherits ``P``'s masks
    when either:

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
    # Resolve each target clade to its founder node; some targets (recombinants,
    # very recent lineages) are not monophyletic and have no founder -- skip them.
    ns_founder = {}
    skipped = []
    for full in ns_full_values:
        node = clade_founders.get(full)
        if node is None:
            skipped.append(full)
        else:
            ns_founder[full] = node
    if skipped:
        lib.log(f"  {len(skipped)} Nextstrain clade(s) have no tree founder node "
                f"(not monophyletic, e.g. recombinant) -- left unmasked: "
                f"{sorted(skipped)}", log_file)

    # node -> Nextstrain clade full-name(s) with that founder (usually one).
    nsnode_to_full = defaultdict(list)
    for full, node in ns_founder.items():
        nsnode_to_full[node].append(full)
    ns_nodes = set(ns_founder.values())

    targets = set(pango_node.values()) | ns_nodes
    if chains is None:
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
            lib.log(f"  note: PANGO clade {p} has no enclosing Nextstrain clade in "
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
        # Union the site sets and mutation lists of every assigned PANGO clade.
        sites = set()
        mutations = set()
        for p in pango_clades:
            sites.update(spec[p]["sites"])
            mutations.update(spec[p]["mutations"])
        clade_masks[ns_clade] = {
            "sites": sites,
            "mutations": mutations,
            "sources": sorted(pango_clades),  # provenance, emitted as a comment
        }
        lib.log(f"  {ns_clade}: {len(sites)} sites, {len(mutations)} mutations "
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
        # m[1:-1] is the position of a refPOSalt mutation string, e.g. "670" of "G670T".
        in_muts = {m for m in block["mutations"] if start <= int(m[1:-1]) <= end}
        dropped_muts += len(block["mutations"]) - len(in_muts)
        block["mutations"] = in_muts
    lib.log(f"Restricted to site_include_range [{start}, {end}]: dropped "
            f"{dropped_sites} UTR/out-of-range site entries and {dropped_muts} "
            f"mutation entries.", log_file)
    return dropped_sites, dropped_muts


def build_masks_for_nodes(spec, node_map, target_full_values, clade_founders,
                          sample_paths_file, start, end, log_file=None, chains=None):
    """Run the full over-mask assignment for one anchoring, for one target level.

    Convenience wrapper around ``assign_masks_by_founder`` -> ``build_clade_masks``
    -> ``restrict_to_include_range``. ``node_map`` maps ``pango_clade -> anchor node
    id``; ``target_full_values`` are the target clade labels (Nextstrain full values,
    or the PANGO clade names for the PANGO-level output). Pass ``chains`` to reuse a
    single sample-paths scan. Returns the in-range
    ``{bare_target_clade: {"sites", "mutations", "sources"}}`` map.
    """
    ns_to_pango = assign_masks_by_founder(
        spec, node_map, target_full_values, clade_founders, sample_paths_file,
        log_file, chains=chains)
    clade_masks = build_clade_masks(spec, ns_to_pango, log_file)
    restrict_to_include_range(clade_masks, start, end, log_file)
    return clade_masks


# ---------------------------------------------------------------------------
# Anchor cross-check: representative vs. founder patristic distance
# ---------------------------------------------------------------------------


def compute_anchor_distances(founder_node, rep_node, sample_paths_file,
                             log_file=None, chains=None):
    """Patristic distance between each clade's founder anchor and representative.

    For every clade in ``founder_node`` returns an entry describing how far its
    founder (MRCA) anchor sits from its representative anchor, walking the tree:

    * ``rep_absent`` -- the representative is not in the tree (no distance).
    * otherwise ``relation`` (identical / one ancestral to the other / siblings via
      a common ancestor), ``d_nodes`` and ``d_muts`` (patristic distances = edges /
      nucleotide mutations summed along both branches from their last common
      ancestor), and ``diff_muts`` (the sorted set of tree mutations separating the
      two anchors).

    Pass precomputed ``chains`` (from ``capture_token_chains``, covering every
    founder + representative node) to reuse a single sample-paths scan; otherwise
    one mutation-preserving scan is done here covering all anchors.
    """
    if chains is None:
        targets = set(founder_node.values()) | {n for n in rep_node.values() if n}
        chains = capture_token_chains(sample_paths_file, targets, log_file)

    def ids_of(chain):
        return [t.split(":", 1)[0] for t in chain]

    out = {}
    for clade, fo in founder_node.items():
        r = rep_node.get(clade)
        if r is None:
            out[clade] = {"founder": fo, "rep": None, "rep_absent": True}
            continue
        rc, fc = chains.get(r), chains.get(fo)
        if not rc or not fc:
            out[clade] = {"founder": fo, "rep": r, "rep_absent": False,
                          "error": "anchor node not found in any sample path"}
            continue
        # Walk both root->anchor chains forward to the last common ancestor (index k);
        # what remains on each side is the branch from that ancestor to each anchor.
        rids, fids = ids_of(rc), ids_of(fc)
        k = 0
        while k < len(rids) and k < len(fids) and rids[k] == fids[k]:
            k += 1
        seg_r, seg_f = rc[k:], fc[k:]  # branches below the last common ancestor
        diff = []
        for t in seg_r:
            diff += _token_muts(t)
        for t in seg_f:
            diff += _token_muts(t)
        if not seg_r and not seg_f:
            relation = "identical"
        elif not seg_f:
            relation = "founder ancestral to representative"
        elif not seg_r:
            relation = "representative ancestral to founder"
        else:
            relation = "siblings (via common ancestor)"
        out[clade] = {
            "founder": fo, "rep": r, "rep_absent": False, "relation": relation,
            # d_muts sums mutations along both branches (len, not set) = patristic;
            # diff_muts dedups them only for display.
            "d_nodes": len(seg_r) + len(seg_f), "d_muts": len(diff),
            "diff_muts": sorted(set(diff)),
        }
    return out


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


# ---------------------------------------------------------------------------
# Emit a clade-keyed masking YAML
# ---------------------------------------------------------------------------


def write_output(out_path, clade_masks, spec, provenance, dropped_exclusions,
                 level="Nextstrain"):
    """Write a clade-keyed masking YAML, one block per clade.

    ``level`` is ``"Nextstrain"`` (the mat-pcp output, over-mask-coarsened to
    Nextstrain clades) or ``"PANGO"`` (the per-PANGO-lineage companion output).
    """
    coarsen = ("and coarsened to Nextstrain clades (over-mask policy for PANGO "
               "clades finer than a Nextstrain clade)."
               if level == "Nextstrain" else
               "keyed by PANGO lineage (each clade accumulates its ancestral PANGO "
               "clades' masks; no Nextstrain coarsening).")
    lines = []
    lines.append("# ============================================================"
                 "=================")
    lines.append("# GENERATED FILE -- DO NOT EDIT BY HAND.")
    lines.append("#")
    lines.append(f"# {level}-clade masking file produced by the s2-clade-mask generator")
    lines.append("# from the PANGO-lineage-keyed UCSC spec data/branchSpecificMask.yml, resolved")
    lines.append("# against a full gisaid-and-public MAT + metadata. Masks are accumulated down the")
    lines.append(f"# tree (topology-based inheritance) {coarsen}")
    lines.append("# Sites and reversion mutations outside the pipeline's")
    lines.append("# site_include_range (the UTRs) are omitted for simplicity -- they")
    lines.append("# are discarded by the range filter at runtime regardless.")
    lines.append("#")
    for line in provenance:
        lines.append(f"# {line}")
    if dropped_exclusions:
        lines.append("#")
        lines.append("# NOTE: the following per-sample 'exclusions' from the source spec")
        lines.append("# cannot be represented at clade granularity and were")
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
        if not sites and not mutations:  # nothing in range for this clade
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

    # CLI args win over config defaults for every path this rule reads/writes.
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
    ns_full_values = load_nextstrain_clade_values(
        cfg["metadata"], cfg["clade_column"], log_file)
    start, end = cfg["site_include_start"], cfg["site_include_end"]

    # The output is founder-anchored; resolve_founder already hard-errors on any
    # clade missing a founder, so nodes_founder resolves every spec clade.
    founder_node, _, _ = lib.read_nodes(nodes_founder)
    rep_node, _, rep_unresolved = lib.read_nodes(nodes_rep)
    lib.log(f"Founder anchors: {len(founder_node)}; representative anchors: "
            f"{len(rep_node)} ({len(rep_unresolved)} absent)", log_file)

    # --- Single sample-paths scan feeding both the anchor cross-check and the
    #     over-mask assignment. Token chains (node:mutations) subsume the id-only
    #     ancestor chains the assignment needs, so one pass over the ~19 GB file
    #     covers every anchor (founder + representative) and Nextstrain founder node.
    ns_founder_nodes = {clade_founders[v] for v in ns_full_values
                        if v in clade_founders}
    all_targets = (set(founder_node.values())
                   | {n for n in rep_node.values() if n}
                   | ns_founder_nodes)
    lib.log("\n=== Scanning sample paths (shared: cross-check + both outputs) ===",
            log_file)
    token_chains = capture_token_chains(sample_paths, all_targets, log_file)
    chains = {n: [t.split(":", 1)[0] for t in toks]
              for n, toks in token_chains.items()}

    # --- Anchor cross-check table (report-only here; enforced in check_anchors) --
    lib.log("\n=== Cross-checking founder vs representative anchors ===", log_file)
    distances = compute_anchor_distances(founder_node, rep_node, sample_paths,
                                         log_file, chains=token_chains)
    write_distance_report(out_comparison, distances, cfg["max_anchor_distance"],
                          set(cfg["representative_absent_ok"]))
    write_distance_data(out_distances, distances)
    lib.log(f"Wrote anchor distance table to {out_comparison} / {out_distances}",
            log_file)

    lib.log("\n=== Building Nextstrain-clade masks (founder-anchored) ===", log_file)
    ns_masks = build_masks_for_nodes(
        spec, founder_node, ns_full_values, clade_founders, sample_paths,
        start, end, log_file, chains=chains)

    lib.log("\n=== Building PANGO-clade masks (founder-anchored) ===", log_file)
    # For the PANGO output every clade is its own target, so rule (b) never
    # coarsens -- each clade just accumulates its ancestral spec clades via rule (a).
    pango_masks = build_masks_for_nodes(
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
    write_output(out_yaml, ns_masks, spec, base_prov, dropped_exclusions,
                 level="Nextstrain")
    lib.log(f"Wrote {out_yaml}", log_file)
    write_output(out_pango, pango_masks, spec, base_prov, dropped_exclusions,
                 level="PANGO")
    lib.log(f"Wrote {out_pango}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
