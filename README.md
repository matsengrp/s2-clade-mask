# s2-clade-mask

Generate `usher_masked_sites.yaml` — the per-**Nextstrain-clade** site/mutation
masking file used by the [`mat-pcp`](https://github.com/matsengrp/mat-pcp)
SARS-CoV-2 pipeline — from the UCSC PANGO-lineage-keyed branch-specific masking
spec (`branchSpecificMask.yml`).

This repository is **standalone**: its one job is to produce
`results/usher_masked_sites.yaml`. Generation needs the full ~8M-sample
gisaid-and-public mutation-annotated tree (MAT) and `matUtils`; application does
not. When the masks are refreshed here, the resulting YAML is copied into
`mat-pcp` (`data/usher_masked_sites.yaml`) — that copy is the only thing `mat-pcp`
needs from this repo.

## Contents

```
s2-clade-mask/
├── README.md
├── CLAUDE.md                       # guidance for Claude Code / contributors
├── config.yaml                     # all run parameters (paths, options)
├── environment.yml                 # conda env (usher/matUtils, python, pyyaml)
├── scripts/
│   └── generate_usher_masked_sites.py   # the generator
├── data/
│   └── branchSpecificMask.yml      # vendored UCSC input spec (pinned commit)
└── results/
    └── usher_masked_sites.yaml     # generated output (tracked)
```

## The problem the generator solves

Two masking systems with incompatible keys:

- **mat-pcp (output):** `usher_masked_sites.yaml`, keyed by **Nextstrain clade**
  (`21A`, `22B`). Masking is by exact string match — no tree involved at runtime.
  Inheritance is *baked in* by repeating the accumulated sites across every
  descendant clade's block.
- **UCSC (input) `branchSpecificMask.yml`:** keyed by **PANGO lineage** (`BA.4`,
  `BQ.1.1`), each anchored to a **representative sample**. Upstream,
  `matUtils mask` applies each mask to that sample's node *and its entire
  subtree*, so inheritance is implicit in tree topology. This spec replaced the
  old `maskDelta.sh` approach; see
  [yatisht/usher#324](https://github.com/yatisht/usher/issues/324).

`generate_usher_masked_sites.py` is an **offline** bridge: it reads the PANGO
spec, uses a full **gisaid-and-public** MAT to compute real tree ancestry (the
GISAID tree is required because several representative samples are GISAID-only),
coarsens PANGO → Nextstrain, and emits the flat Nextstrain-keyed file `mat-pcp`
already knows how to consume.

## Installation

```bash
conda env create -f environment.yml
conda activate s2-clade-mask
# sanity check
python scripts/generate_usher_masked_sites.py --help
```

`usher` supplies the `matUtils` binary; the generator itself only needs the
Python standard library plus `pyyaml`.

## Usage

All run parameters live in `config.yaml`. Edit it to point at your input tree /
metadata, then run the generator **from the repo root** (so relative paths
resolve):

```bash
python scripts/generate_usher_masked_sites.py --config config.yaml
```

Use the **full gisaid-and-public** MAT — not the public-only tree (it is missing
several GISAID-only representatives) and not a filtered per-run tree (all clades
must be present / annotated). At the Hutch a recent tree lives at
`/fh/fast/matsen_e/shared/viral-dasm/mat-trees/<date>-gisaid-and-public/`; point
`config.yaml` at the **GISAID-named** files (`*.gisaidNames.pb.gz` /
`*.gisaidNames.metadata.tsv.gz`), whose sample names are what the spec
representatives match against.

`config.yaml` keys:

| Key | Meaning |
|---|---|
| `tree`, `metadata` | full gisaid-and-public MAT `.pb.gz` and its metadata (required unless `sample_paths` + `clade_paths` are precomputed) |
| `branch_spec` | vendored UCSC spec (`data/branchSpecificMask.yml`) |
| `output` | where to write the YAML (`results/usher_masked_sites.yaml`) |
| `old_yaml` | reference file for the superset-recovery check (optional) |
| `clade_column` | metadata clade column (default `Nextstrain_clade`) |
| `site_include_start` / `site_include_end` | coding range `[266, 29674]`; UTR masks outside it are dropped |
| `strict` | fail if old masks are not fully recovered |
| `union_old` | force-include any unrecovered old masks |
| `matutils`, `workdir` | binary name and scratch dir for the large paths files |
| `sample_paths`, `clade_paths` | precomputed matUtils outputs (skip the extract passes) |
| `log` | optional run-log file |
| `representative_name_overrides` | spec rep name → exact tree name, for reviewed renames (see §3) |
| `representative_fallbacks` | clade → approved fallback when its rep is absent (see §3) |

The `config.yaml` also carries the human-approved `representative_name_overrides`
/ `representative_fallbacks` tables (see §3). Because the config is
version-controlled, those entries are reviewed exactly like code.

`matUtils extract --sample-paths` on an ~8M-sample tree produces a large paths
file (tens of GB); the generator reads it only until every anchor node's ancestor
chain has been seen (early exit), so memory stays modest. Point `workdir` at a
disk with room, or set `sample_paths` / `clade_paths` to skip `matUtils`.

> **Note — GISAID redistribution.** Outputs derived from a GISAID tree inherit
> GISAID's redistribution terms. `usher_masked_sites.yaml` itself contains only
> site numbers and mutation strings (no sequence data), but representative sample
> names appear in the run log.

## How the script works, stage by stage

### 1. Parse the UCSC spec (`parse_branch_spec`)

For each of the 25 PANGO clades it collects the `representative`
(+ `representativeBacktrack`), expands every `ranges: [[lo, hi]]` into individual
integer positions and merges with `sites` → a flat **site set**; validates each
`reversions` string against `^[ACGT][0-9]+[ACGT]$` → a **mutation list**; keeps
`exclusions` aside. It also asserts the spec has a `version` key (fails loudly
if the upstream format changes).

### 2. Get node IDs from the tree

Two `matUtils extract` passes on the full `.pb.gz`:

- `--sample-paths` → each sample's root → tip node chain. This is the **primary**
  node source: each PANGO clade is anchored at its *representative sample's* node
  (`resolve_representative_nodes`), the same node `matUtils mask` applies the mask
  to upstream. It is also read to obtain each anchor node's ancestor chain.
- `--clade-paths` → `{clade: founder_node_id}` — the founder (MRCA) node of
  *every* annotated clade, both PANGO (`BA.4`) and Nextstrain (`"22A (BA.4)"`)
  (`parse_clade_paths`). Nextstrain founders anchor the Nextstrain side of the
  over-mask rule (§5); PANGO founders are used only as an approved per-clade
  fallback (§3).

### 3. Resolve each PANGO clade to a node (`resolve_clade_nodes`)

Each PANGO clade is anchored at its **representative sample's node**, resolved
from `--sample-paths` with the exact upstream `branchSpecificMask.py` logic (take
the last path node, back up past the sample's own leaf, then back up
`representativeBacktrack` more). This is faithful to how `matUtils mask` applies
each mask — to the representative's node *and its subtree* — rather than to the
clade's MRCA, which can sit higher and over-broaden the mask.

Because several representatives are **GISAID-only** names absent from the
public-only tree, the generator must run against the full **gisaid-and-public**
MAT. Against the 2026-04-22 tree, 22 of the 25 PANGO clade representatives match
by exact name, two need an approved name override (`B.1.617.2`, `XBB`), and one
(`P.1`) is absent and uses an approved fallback — so all 25 resolve.

Resolution is deliberately **explicit and auditable** — a person signs off on
every deviation:

1. **Representative sample (primary).** Match `representative` **exactly** against
   the tree sample name *before the first `|`*. Never a fuzzy or substring match:
   e.g. `Foo/xGBRC714bx/2021` must not resolve `IND/GBRC714b/2021`.
2. **Approved name override** (config `representative_name_overrides`). When the
   spec name does not match the tree verbatim, a maintainer records the confirmed
   tree name. Currently two entries, both pure country-abbreviation differences
   verified against the 2026-04-22 tree: `IND/GBRC714b/2021` →
   `India/GJ-GBRC714b/2021` (clade `B.1.617.2`) and `BGD/icddrb_TND_06_1053/2022`
   → `Bangladesh/icddrb_TND_06_1053/2022` (clade `XBB`).
3. **Approved fallback** (config `representative_fallbacks`). When a representative
   is genuinely absent, a maintainer approves an alternative per clade —
   `via: founder` (the annotated MRCA node) or `via: sample` (an approved alternate
   sample). Currently one entry: `P.1`, whose representative
   `FRA/IHUCOVID-005193-N1/2021` is not in the 2026-04-22 tree, falls back to its
   `20J (Gamma)` founder. (P.1's only mask is the deletion `[11288, 11296]`, also
   masked via `B.1.1.7` and `B.1.351`, so nothing is lost.)
4. **Otherwise, hard error.** A representative that is neither found nor covered by
   an approved override/fallback stops the run, naming the clade and instructing
   the maintainer to add an entry — so a future spec revision that renames or drops
   a representative can never be silently mis-resolved.

Both tables live in `config.yaml` (`representative_name_overrides` /
`representative_fallbacks`), validated at load time; revisit them whenever the
spec or the tree is refreshed.

### 4. Enumerate Nextstrain clades (`load_nextstrain_clade_values`)

From the metadata `Nextstrain_clade` column. Keeps full values (`"22A (BA.4)"`)
for founder lookup; the **output key** is the bare token (`"22A"`) via
`bare_clade`, exactly matching mat-pcp's `get_clades` awk split.

### 5. Assign PANGO masks to Nextstrain clades (`assign_masks_by_founder`)

The core logic. Each PANGO clade `P` is anchored at its representative node
`r(P)` (§3); each Nextstrain clade `C` at its founder (MRCA) node. Working from
these single anchor nodes (not individual member samples) is deliberate: it is
noise-free, whereas a per-sample membership rule lets recombinants and mislabeled
samples bleed masks across unrelated clades. A Nextstrain clade `C` inherits PANGO
clade `P`'s masks iff:

- **(a)** `r(P)` lies on the root → `founder(C)` chain — `P` is an ancestor
  of `C`, so all of `C` descends from `r(P)`; the mask applies **exactly**; or
- **(b)** `C` is the *deepest* Nextstrain clade enclosing `r(P)` (last
  Nextstrain founder on the root → `r(P)` chain) — `P` is finer than any
  Nextstrain clade, so its mask is **widened** to all of `C` (over-mask
  coarsening).

Rule (b) attributes a fine PANGO clade to the clade its samples *actually carry* —
e.g. `KP.3.1.1` → `24E`, not its ancestor `24A`; `21L` (BA.2) gets only BA.2's own
masks. Clades with no founder node (non-monophyletic `recombinant`) are logged and
left unmasked.

### 6. Accumulate & emit (`build_clade_masks` → `write_output`)

Union the site sets and mutation lists of all PANGO clades assigned to each
Nextstrain clade; write one `mask_<bare_clade>` block with `clades` / `sites` /
`mutations`, a `# accumulated from PANGO clades: [...]` provenance comment per
block, and a file header (source commit, tree, metadata, dropped-exclusion note).

### 7. Drop UTR / out-of-range masks (`restrict_to_include_range`)

Sites and reversion mutations outside `site_include_range` (`[266, 29674]` by
default — everything outside is the 5′/3′ UTR) are removed. `mat-pcp` already
discards those positions via its range filter at runtime, so masking them was
redundant; omitting them keeps the file focused on the masks that actually affect
results. The old-scheme recovery check (below) treats old UTR masks as
intentionally dropped rather than as failures.

## Output schema

`usher_masked_sites.yaml` is a flat map of `mask_<clade>` blocks:

```yaml
mask_22B:
  clades: [22B]                 # Nextstrain clade(s) this block applies to
  sites: [11288, 11289, ...]    # whole-site masks (1-based genome positions)
  mutations: ['G670T', ...]     # optional; directional refPOSalt exclusions
```

Mapping of source fields:

- `ranges` + `sites` → **site masks**.
- `reversions` → **mutation masks** (directional `refPOSalt`, matching mat-pcp's
  `nt_mutation`, *not* expanded to whole sites).
- `exclusions` → dropped + logged (cannot express per-sample un-masking at
  clade granularity).

**Consumed in `mat-pcp`** by `count_mutations.py` (unions `mutations` into its
`muts_to_exclude` set) and `aggregate_mutations_to_exclude.py` (emits one row per
mutation). Blocks without a `mutations:` field behave exactly as before, so the
field is backward-compatible.

## How the new masks match the old ones

> **Note — figures pending refresh.** The counts below were measured against a
> **public-only** tree (`public-2026-07-01`) under an earlier *founder-primary*
> resolution. The current *representative-primary* logic (§3) and the switch to a
> gisaid-and-public tree will shift the exact numbers on the next regeneration;
> the acceptance criterion (old-scheme superset recovery) still holds. Re-measure
> when `results/usher_masked_sites.yaml` is regenerated.

The generator enforces an explicit acceptance criterion (`check_recovers_old`):
the new file must be a **superset** — for every `(clade, site)` in the old file
*within the `site_include_range`*, that site must still be masked for that clade
(as a whole-site mask, or now at the **mutation** level — i.e. ≥ 1 new mutation at
that position for that clade, an intentional, more precise "downgrade"). Old masks
in the UTRs are intentionally not recovered (§7) and are reported separately.

With the recovery check enabled (`old_yaml` set) on the public-2026-07-01 snapshot,
the run reported full recovery (exit 0): every old `(clade, site)` pair preserved,
141 whole-site → mutation downgrades (in `22B` / `22D` / `22E`), zero truly
missing, plus large new coverage for post-2022 clades the hand-made file never had
(14 → 43 clades). `strict: true` makes any *true* shortfall a hard error;
`union_old: true` would force-include unrecovered old masks.

## Updating the vendored spec

Re-download `branchSpecificMask.yml` from the UCSC kent repo
(`src/hg/utils/otto/sarscov2phylo/branchSpecificMask.yml`) at a newer commit,
update the provenance block at the top of `data/branchSpecificMask.yml`
(current: commit `2c7413b`, 2026-03-13), then regenerate as above. If the run
fails with an unresolved-representative error, a representative was renamed or
dropped: confirm the sample by hand in the tree and add an approved entry to
`representative_name_overrides` (renamed) or `representative_fallbacks` (absent)
in `config.yaml` — the script never guesses.

## Caveats / accepted approximations

- **Coarsening loss (accepted).** PANGO clades finer than any Nextstrain clade
  are widened to the enclosing Nextstrain clade (over-mask). Logged per clade.
- **Reversions.** Handled at the mutation level, not blanket-masked over the site.
- **Exclusions.** Per-sample un-masking of a site cannot be expressed at
  Nextstrain granularity → skipped + logged (currently only BA.2.86 site 23599).
- **Recombinants.** The non-monophyletic `recombinant` clade has no single
  founder node and is left unmasked.
- **Representative resolution (manual sign-off).** Representatives match the tree
  by exact pre-`|` name only; overrides and fallbacks each require a hard-coded,
  human-reviewed entry (see §3).
- **GISAID tree required.** Several representatives are GISAID-only; outputs
  inherit GISAID redistribution terms (the YAML itself carries no sequence data).
- **Drift.** The UCSC spec grows as new lineages are designated; the vendored copy
  is pinned to a specific commit and should be refreshed periodically.
