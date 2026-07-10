#!/usr/bin/env python3
"""Pipeline rule entrypoint: resolve PANGO clades to tree nodes in one anchoring mode.

Two modes, written as separate node-map artifacts (``s2clademask.write_nodes``
schema) so the downstream ``compare_and_build`` rule can diff them and build the
final output:

* ``--mode representative`` (Mode 1, string match): anchor each PANGO clade at its
  representative sample's node via an exact pre-``|`` name match (applying the
  approved ``representative_name_overrides``, but **no** founder fallback). Clades
  whose representative is absent from the tree are recorded as ``unresolved``.
* ``--mode founder`` (Mode 2): anchor each PANGO clade at its annotated founder
  (MRCA) node from ``--clade-paths``. Clades with no annotated founder are a hard
  error (the founder is the authoritative anchor and must exist).

The resolution logic itself lives in this script (it is used by no other rule); the
shared parsers / config / node-map schema come from ``s2clademask``.

    python scripts/resolve_nodes.py -c config.yaml --mode representative \\
        --sample-paths results/sample_paths.tsv --out results/nodes_representative.yaml
    python scripts/resolve_nodes.py -c config.yaml --mode founder \\
        --clade-paths results/clade_paths.tsv --out results/nodes_founder.yaml
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2clademask as lib


# ---------------------------------------------------------------------------
# Mode 1: representative-sample resolution (string match)
# ---------------------------------------------------------------------------


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
    remaining = set(wanted)  # names still to find; drives the early exit below
    with lib.open_maybe_gzip(sample_paths_file) as f:
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
                # If the last node is the sample's own leaf, step past it first.
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
            if not remaining:  # every wanted sample seen -- stop scanning the huge file
                break

    return clade_to_node


def resolve_representative_only(spec, sample_paths_file, name_overrides,
                                log_file=None):
    """Mode 1 (string match): anchor each clade at its representative sample node.

    Steps 1-2 of the resolution policy, with **no** founder/sample fallback:

    1. **Representative sample.** Resolve ``spec[clade]["representative"]``
       from the sample-paths file, matching the name *exactly* before the first
       ``|`` (upstream ``branchSpecificMask.py`` logic, respecting
       ``representativeBacktrack``). This anchors the mask at the same node
       ``matUtils mask`` would.
    2. **Approved name override.** If the spec name does not match the tree
       verbatim but a human has recorded the correct tree name in
       ``representative_name_overrides``, resolve using that name.

    Returns ``(clade_to_node, via, unresolved)``: the resolved nodes, how each was
    resolved (``"representative"`` / ``"name override"``), and the sorted list of
    clades whose representative is absent from the tree (each of which must be
    listed in ``representative_absent_ok`` or the run fails). Run against the full
    **gisaid-and-public** MAT; the public-only tree is missing several GISAID-only
    representatives.
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
        lib.log(f"  name override (approved): {clade} '{spec_name}' -> '{tree_name}'",
                log_file)
    return clade_to_node, via, unresolved


# ---------------------------------------------------------------------------
# Mode 2: founder (MRCA) resolution
# ---------------------------------------------------------------------------


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


def resolve_founder_or_fail(spec, clade_founders, log_file=None):
    """Founder anchor for every PANGO clade; hard error if any is missing.

    Returns ``clade_to_node`` mapping each spec clade to its annotated founder
    (MRCA) node. A clade with no annotated founder in the tree is a hard error
    (design invariant #1: the founder is the anchor and must be present -- run
    against the full gisaid-and-public MAT, where every spec clade is annotated).
    """
    founder_node, no_founder = resolve_founder_nodes(spec, clade_founders)
    if no_founder:
        detail = "\n".join(
            f"    {c}: no annotated founder node in the tree" for c in no_founder)
        raise ValueError(
            "The following PANGO clade(s) have no annotated founder (MRCA) node in "
            "the tree, so they cannot be founder-anchored. Run against the full "
            "gisaid-and-public MAT where every clade is annotated:\n" + detail)
    for clade in sorted(founder_node):
        lib.log(f"  {clade} -> {founder_node[clade]} (founder)", log_file)
    return founder_node


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", "-c", default="config.yaml",
                   help="Path to the YAML config (default: config.yaml).")
    p.add_argument("--mode", choices=["representative", "founder"], required=True,
                   help="Anchoring mode to resolve.")
    p.add_argument("--sample-paths", help="matUtils --sample-paths TSV (mode representative).")
    p.add_argument("--clade-paths", help="matUtils --clade-paths TSV (mode founder).")
    p.add_argument("--out", required=True, help="Output node-map YAML path.")
    cli = p.parse_args(argv)

    # Load + validate the config, then parse the spec both modes anchor from.
    cfg = lib.load_config(cli.config)
    log_file = cfg["log"]
    lib.validate_resolution_overrides(cfg["representative_name_overrides"],
                                      cfg["representative_absent_ok"])
    spec, version = lib.parse_branch_spec(cfg["branch_spec"])
    lib.log(f"Parsed {len(spec)} PANGO clades, spec version {version}", log_file)

    # Make sure the output directory exists before we write the node map.
    out_dir = os.path.dirname(os.path.abspath(cli.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Each mode produces the same node-map artifact from a different input file and
    # anchoring rule; `via` records how each clade was anchored, `meta` is provenance.
    if cli.mode == "representative":
        # Mode 1 reads the (large) per-sample paths and string-matches representatives.
        sample_paths = cli.sample_paths
        if not sample_paths:
            raise ValueError("mode representative needs --sample-paths")
        lib.log("Mode 1 (representative / string match): resolving representative "
                "sample nodes (for the founder cross-check)...", log_file)
        clade_to_node, via, unresolved = resolve_representative_only(
            spec, sample_paths, cfg["representative_name_overrides"], log_file)
        # A representative missing from the tree is only acceptable if the maintainer
        # has approved it in representative_absent_ok (its cross-check is then
        # skipped). Any other absence is a hard error naming the clade -- never a
        # silent guess -- so a spec/tree rename can't slip through unnoticed.
        not_ok = [c for c in unresolved if c not in cfg["representative_absent_ok"]]
        if not_ok:
            detail = "\n".join(
                f"    {c}: representative '{spec[c]['representative']}' not found"
                for c in sorted(not_ok))
            raise ValueError(
                "The following clade(s) have a representative absent from the tree "
                "and are not listed in config representative_absent_ok. Confirm the "
                "sample by hand, then add a representative_name_overrides entry (if "
                "renamed) or list the clade in representative_absent_ok (if truly "
                f"absent):\n{detail}")
        meta = {"mode": "representative", "spec_version": version,
                "source": os.path.abspath(sample_paths)}
    else:
        # Mode 2 reads the small clade-paths TSV and takes each clade's founder node.
        clade_paths = cli.clade_paths
        if not clade_paths:
            raise ValueError("mode founder needs --clade-paths")
        clade_founders = lib.parse_clade_paths(clade_paths)
        lib.log(f"Parsed {len(clade_founders):,} clade annotations from {clade_paths}",
                log_file)
        lib.log("Mode 2 (founder / MRCA): resolving annotated founder nodes "
                "(each clade's anchor)...", log_file)
        # The founder is the authoritative anchor, so it must exist for every clade;
        # resolve_founder_or_fail raises (no unresolved) if any is missing.
        clade_to_node = resolve_founder_or_fail(spec, clade_founders, log_file)
        via = {c: "founder" for c in clade_to_node}
        unresolved = []
        meta = {"mode": "founder", "spec_version": version,
                "source": os.path.abspath(clade_paths)}

    # Emit the node map (shared schema) for the downstream compare_and_build rule.
    lib.write_nodes(cli.out, clade_to_node, via=via, unresolved=unresolved, meta=meta)
    lib.log(f"Resolved {len(clade_to_node)}/{len(spec)} clades ({len(unresolved)} "
            f"unresolved) -> {cli.out}", log_file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
