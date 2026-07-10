# s2-clade-mask

`s2-clade-mask` generates two SARS-CoV-2 site/mutation masking files from the
UCSC PANGO-lineage-keyed branch-specific masking spec (`branchSpecificMask.yml`):

- **`results/usher_masked_sites.yaml`** — keyed by **Nextstrain clade** (`21A`,
  `22B`). Consumed by the [`mat-pcp`](https://github.com/matsengrp/mat-pcp)
  SARS-CoV-2 pipeline.
- **`results/usher_masked_sites_pango.yaml`** — keyed by **PANGO lineage** (`BA.4`,
  `BQ.1.1`). Used in other contexts that work at PANGO-lineage granularity.

Both carry the same underlying masks and are produced together from a single run;
they differ only in the clade level they are keyed by — the Nextstrain file
coarsens PANGO clades finer than a Nextstrain clade up to the enclosing Nextstrain
clade (over-mask), while the PANGO file keeps them at full lineage resolution.

This repository is **standalone**. Generation needs the full gisaid-and-public
mutation-annotated tree (MAT) — ~17M samples on the 2026-04-22 tree, and growing —
plus `matUtils`; downstream application of either file does not. When the masks are
refreshed here, the Nextstrain-clade YAML is copied into `mat-pcp`
(`data/usher_masked_sites.yaml`) — that copy is the only thing `mat-pcp` needs from
this repo.

## Contents

```
s2-clade-mask/
├── README.md
├── CLAUDE.md                       # guidance for Claude Code / contributors
├── config.yaml                     # all run parameters (paths, options)
├── environment.yml                 # conda env (usher/matUtils, snakemake, python, pyyaml, gh)
├── Snakefile                       # staged pipeline (primary run path)
├── scripts/
│   ├── s2clademask.py                   # shared library (all generation logic)
│   ├── extract_paths.py                 # rule: one matUtils extract pass
│   ├── resolve_nodes.py                 # rule: resolve clades in one anchoring mode
│   ├── compare_and_build.py             # rule: emit both YAMLs + anchor cross-check table
│   ├── check_anchors.py                 # rule: enforce max_anchor_distance
│   └── check_recovery.py                # rule: superset recovery vs old scheme
├── data/
│   ├── branchSpecificMask.yml           # vendored UCSC input spec (pinned commit)
│   └── usher_masked_sites_old_scheme.yml # vendored jbloomlab old-scheme baseline (pinned commit)
└── results/
    ├── usher_masked_sites.yaml       # generated Nextstrain-clade output (tracked)
    └── usher_masked_sites_pango.yaml # generated PANGO-clade output (tracked)
```

## The problem the generator solves

The UCSC spec `branchSpecificMask.yml` describes SARS-CoV-2 masking as clade-by-clade
*logic* over the UShER tree, keyed by PANGO lineage: for each clade it lists the
sites and reversion mutations to mask, anchored at that clade's node in the tree.
Upstream, `matUtils mask` applies each clade's masks to that node **and its entire
subtree**, so a clade's true mask is its own entries plus everything inherited from
its ancestors — that inheritance lives in the tree topology and is never written
out. (This spec replaced the old `maskDelta.sh` approach; see
[yatisht/usher#324](https://github.com/yatisht/usher/issues/324).)

The goal of this pipeline is to **parse that logic into an explicit, complete list
of masked sites and mutations for each clade.** It resolves every clade to its tree
node, walks the topology to accumulate each clade's own masks together with all the
masks it inherits, and writes the result out as the two flat, self-contained YAMLs
above — so a consumer can mask by clade name alone, with no tree needed at runtime.

## Installation

```bash
conda env create -f environment.yml
conda activate s2-clade-mask
snakemake -n   # sanity check: dry-run the pipeline DAG
```

The `s2-clade-mask` conda env created above supplies everything the pipeline
needs: `usher` (the `matUtils` binary), `snakemake` (which runs the staged
pipeline), `gh` (GitHub CLI, for issues/PRs), and Python with `pyyaml`. The
generation code itself only uses the Python standard library plus `pyyaml`.

## Usage

All run parameters live in `config.yaml`. Edit it to point at your input tree /
metadata, then run **from the repo root** (so the relative paths resolve):

```bash
snakemake --cores 4
```

The generation is decomposed into cached [Snakemake](https://snakemake.readthedocs.io)
rules (see [How the Snakemake pipeline works](#how-the-snakemake-pipeline-works)),
so each expensive stage re-runs only when its inputs change — the two
`matUtils extract` passes especially (the `--sample-paths` pass is ~1 h and produces
a ~19 GB file).

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
| `tree`, `metadata` | full gisaid-and-public MAT `.pb.gz` and its metadata (both required) |
| `branch_spec` | vendored UCSC spec (`data/branchSpecificMask.yml`) |
| `output` | Nextstrain-clade YAML (`results/usher_masked_sites.yaml`) |
| `output_pango` | PANGO-clade companion YAML (`results/usher_masked_sites_pango.yaml`) |
| `max_anchor_distance` | fail if a clade's founder anchor is more than this many mutations (patristic) from its representative (default `4`; see §3) |
| `old_yaml` | reference file for the superset-recovery check (optional); the vendored original scheme, `data/usher_masked_sites_old_scheme.yml` |
| `clade_column` | metadata clade column (default `Nextstrain_clade`) |
| `site_include_start` / `site_include_end` | coding range `[266, 29674]`; UTR masks outside it are dropped |
| `strict` | fail if old masks are not fully recovered |
| `union_old` | force-include any unrecovered old masks |
| `matutils`, `workdir` | binary name and scratch dir for the large paths files |
| `nodes_representative`, `nodes_founder` | founder / representative node-map intermediates (default `results/nodes_{founder,representative}.yaml`) |
| `comparison_report` | human anchor-distance table (default `results/anchoring_comparison.txt`) |
| `anchor_distances` | machine-readable anchor distances (default `results/anchor_distances.yaml`) |
| `anchor_check_report` | anchor cross-check pass/fail result (default `results/anchor_check.txt`) |
| `recovery_report` | old-scheme recovery pass/fail (default `results/recovery_report.txt`) |
| `recovery_summary` | per-clade added-mask table (default `results/recovery_summary.txt`) |
| `log` | optional run-log file (default unset) |
| `representative_name_overrides` | spec rep name → exact tree name, for reviewed renames (see §3) |
| `representative_absent_ok` | clades whose representative is genuinely absent from the tree (cross-check skipped; founder anchoring still applies — see §3) |

All the intermediate/report paths above default under `results/`, are gitignored,
and are regenerated on every run; set one only to relocate that artifact.

The `config.yaml` also carries the human-approved `representative_name_overrides`
/ `representative_absent_ok` tables (see §3). Because the config is
version-controlled, those entries are reviewed exactly like code.

`matUtils extract --sample-paths` on a tree this size (~17M samples) produces a large paths
file (tens of GB); the generator reads it only until every anchor node's ancestor
chain has been seen (early exit), so memory stays modest. Point `workdir` at a
disk with room.

> **Note — GISAID redistribution.** Outputs derived from a GISAID tree inherit
> GISAID's redistribution terms. `usher_masked_sites.yaml` itself contains only
> site numbers and mutation strings (no sequence data), but representative sample
> names appear in the run log.

## How the Snakemake pipeline works

The numbered stages below detail the generation logic; where a stage corresponds to
one or more rules, the rule name(s) and the `results/…` file(s) they write are noted
inline.

### 1. Parse the UCSC spec (`parse_branch_spec`)

For each of the 25 PANGO clades it collects the `representative`
(+ `representativeBacktrack`), expands every `ranges: [[lo, hi]]` into individual
integer positions and merges with `sites` → a flat **site set**; validates each
`reversions` string against `^[ACGT][0-9]+[ACGT]$` → a **mutation list**; and sets
aside any `exclusions` — `{site, representative}` entries naming a sub-lineage that
should be *un-masked* at that site (e.g. `BA.2.86` masks site 23599, but its
`XFJ.3.1` sub-lineage has a real reversion there and is excluded). Such
per-sub-lineage carve-outs can't be expressed at clade granularity, so they are
dropped and logged in §6, leaving the site masked clade-wide.

The spec is a single mapping of a top-level `version` key plus one entry per clade,
so the parser removes `version` first (leaving every other top-level key to be read
as a clade) and errors if it is missing — a cheap check that the upstream file still
has that shape. The `version` value itself is only logged, not compared against an
expected number.

### 2. Get node IDs from the tree

Two `matUtils extract` passes on the full `.pb.gz`:

- `--clade-paths` → `{clade: founder_node_id}` — the founder (MRCA) node of
  *every* annotated clade, both PANGO (`BA.4`) and Nextstrain (`"22A (BA.4)"`)
  (`parse_clade_paths`). PANGO founders are the **primary** anchor for each mask
  (§3); Nextstrain founders anchor the Nextstrain side of the over-mask rule (§5).
- `--sample-paths` → each sample's root → tip node chain. Read to resolve each
  clade's *representative sample* node (`resolve_representative_nodes`) for the
  founder cross-check (§3), and to obtain each anchor node's ancestor chain (§5).

*Rules `extract_sample_paths` → `results/sample_paths.tsv` and `extract_clade_paths`
→ `results/clade_paths.tsv`.*

### 3. Anchor each PANGO clade at its founder, cross-checked by its representative (`resolve_founder_or_fail` + `compute_anchor_distances`)

Each PANGO clade is anchored at its annotated **founder (MRCA) node** (from
`--clade-paths`). In a full gisaid-and-public MAT every spec clade is annotated, so
a **missing founder is a hard error**. The founder is the stable anchor; we do
**not** fall back to the representative.

The clade's **representative sample** is resolved (from `--sample-paths`, exact
upstream `branchSpecificMask.py` logic: last path node, back up past the sample's
own leaf, then `representativeBacktrack` more) only to **cross-check** the founder.
The pipeline asserts that the founder and representative anchors are within
`max_anchor_distance` **mutations** of each other (patristic distance through their
common ancestor); a violation is a hard error (rule `check_anchors`). This keeps
the deliberate divergence from upstream `matUtils mask`'s representative-node choice
provably tiny while gaining the founder's stability.

> **Why founder-primary?** The representative and founder are almost always within
> a few mutations (empirical max **4**, on the 2026-04-22 tree). But a representative
> can sit *just across* a Nextstrain-clade boundary from the founder — e.g. `BA.2.75`,
> whose representative attaches to the `21L` (BA.2) stem one step above the `22D`
> (BA.2.75) crown, over-broadening its coarsened mask onto all of `21L`. Anchoring
> at the founder avoids that; the cross-check catches any case where the two are
> *not* close, which would signal a genuinely mis-resolved anchor.

Representative resolution is deliberately **explicit and auditable** — a person
signs off on every deviation:

1. **Exact match.** Match `representative` **exactly** against the tree sample name
   *before the first `|`*. Never a fuzzy or substring match.
2. **Approved name override** (config `representative_name_overrides`). When the
   spec name does not match the tree verbatim, a maintainer records the confirmed
   tree name. Currently two entries, both pure country-abbreviation differences
   verified against the 2026-04-22 tree: `IND/GBRC714b/2021` →
   `India/GJ-GBRC714b/2021` (clade `B.1.617.2`) and `BGD/icddrb_TND_06_1053/2022`
   → `Bangladesh/icddrb_TND_06_1053/2022` (clade `XBB`).
3. **Approved absent** (config `representative_absent_ok`). When a representative is
   genuinely absent from the tree, the maintainer lists the clade here: the founder
   cross-check is skipped for it (founder anchoring still applies). Currently one
   entry: `P.1`, whose representative `FRA/IHUCOVID-005193-N1/2021` is not in the
   2026-04-22 tree; it is founder-anchored at its `20J (Gamma)` MRCA. (P.1's only
   mask, the deletion `[11288, 11296]`, is also masked via `B.1.1.7` / `B.1.351`.)
4. **Otherwise, hard error.** A representative that is absent and not listed stops
   the run, naming the clade — so a future spec/tree change that renames or drops a
   representative can never be silently ignored.

Both tables live in `config.yaml` (`representative_name_overrides` /
`representative_absent_ok`), validated at load time; revisit them whenever the spec
or the tree is refreshed.

*Rules `resolve_founder` → `results/nodes_founder.yaml` and `resolve_representative`
→ `results/nodes_representative.yaml` produce the two anchors; `compare_and_build`
writes the distance tables (`anchoring_comparison.txt`, `anchor_distances.yaml`) and
`check_anchors` → `anchor_check.txt` enforces the bound (see
[Anchor cross-check](#anchor-cross-check-founder-vs-representative)).*

### 4. Enumerate Nextstrain clades (`load_nextstrain_clade_values`)

From the metadata `Nextstrain_clade` column. Keeps full values (`"22A (BA.4)"`)
for founder lookup; the **output key** is the bare token (`"22A"`) via
`bare_clade`, exactly matching mat-pcp's `get_clades` awk split.

### 5. Assign PANGO masks to target clades (`assign_masks_by_founder`)

The core logic, run for **two target levels**: Nextstrain clades (the mat-pcp
output) and PANGO clades (the companion output). Each PANGO clade `P` is anchored
at its founder node `f(P)` (§3); each target clade `C` at its own founder (MRCA)
node. Working from single anchor nodes (not individual member samples) is
deliberate: it is noise-free, whereas a per-sample membership rule lets recombinants
and mislabeled samples bleed masks across unrelated clades. A target clade `C`
inherits PANGO clade `P`'s masks iff:

- **(a)** `f(P)` lies on the root → `founder(C)` chain — `P` is an ancestor
  of `C`, so all of `C` descends from `f(P)`; the mask applies **exactly**; or
- **(b)** `C` is the *deepest* target clade enclosing `f(P)` (last target founder
  on the root → `f(P)` chain) — `P` is finer than any target clade, so its mask is
  **widened** to all of `C` (over-mask coarsening).

For **Nextstrain** targets, rule (b) attributes a fine PANGO clade to the Nextstrain
clade its samples *actually carry* — e.g. `KP.3.1.1` → `24E`, not its ancestor
`24A`. For **PANGO** targets, rule (b) reduces to self-assignment (a clade encloses
its own founder), so each PANGO clade accumulates its own masks plus those of its
ancestral spec clades (e.g. `BQ.1.1` ⊇ `BQ.1` ⊇ `BA.5`) with no coarsening. Clades
with no founder node (non-monophyletic `recombinant`) are logged and left unmasked.

### 6. Accumulate & emit (`build_clade_masks` → `write_output`)

Union the site sets and mutation lists of all PANGO clades assigned to each target
clade; write one `mask_<bare_clade>` block with `clades` / `sites` / `mutations`, a
`# accumulated from PANGO clades: [...]` provenance comment per block, and a file
header (source commit, tree, metadata, dropped-exclusion note). Done once for the
Nextstrain-clade output and once for the PANGO-clade output.

### 7. Drop UTR / out-of-range masks (`restrict_to_include_range`)

Sites and reversion mutations outside `site_include_range` (`[266, 29674]` by
default — everything outside is the 5′/3′ UTR) are removed. `mat-pcp` already
discards those positions via its range filter at runtime, so masking them was
redundant; omitting them keeps the file focused on the masks that actually affect
results. The old-scheme recovery check (below) treats old UTR masks as
intentionally dropped rather than as failures.

*Stages 4–7 all run in the `compare_and_build` rule, which writes
`results/usher_masked_sites.yaml` and `results/usher_masked_sites_pango.yaml`. A
final rule, `check_recovery` → `results/recovery_report.txt` and
`results/recovery_summary.txt`, checks that Nextstrain output against the old scheme
(see [How the new masks match the old ones](#how-the-new-masks-match-the-old-ones)).*

## Anchor cross-check (founder vs. representative)

Because the output is founder-anchored while upstream `matUtils mask` uses the
representative node, the pipeline resolves each clade **both** ways and quantifies
how far apart the two anchors sit. `compare_and_build` writes the table to
`results/anchoring_comparison.txt` (human) and `results/anchor_distances.yaml`
(machine-readable); for every clade it reports the **patristic distance** between
the founder and representative anchors — in tree nodes and in **nucleotide
mutations** — and the exact set of mutations separating them. `check_anchors` then
**fails the run** if any clade's distance exceeds `max_anchor_distance` (default 4,
the empirical maximum), or if a representative is absent without being listed in
`representative_absent_ok`. On the 2026-04-22 tree every clade is within 4
mutations, so the assertion passes. (This replaces the earlier ad-hoc
`compare_anchoring.py` diagnostic.)

## Output schema

Both YAMLs are a flat map of `mask_<clade>` blocks; they differ only in the key
level — `results/usher_masked_sites.yaml` is keyed by **Nextstrain clade** (for
mat-pcp), `results/usher_masked_sites_pango.yaml` by **PANGO lineage**:

```yaml
mask_22B:
  clades: [22B]                 # clade(s) this block applies to
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

The "old scheme" is the original hand-made
[jbloomlab/SARS2-mut-fitness `usher_masked_sites.yaml`](https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/data/usher_masked_sites.yaml),
vendored at a pinned commit as `data/usher_masked_sites_old_scheme.yml` (see its
provenance header). This is a fixed historical baseline the acceptance criterion
checks against — always re-run against this vendored file, never against a
previous `results/usher_masked_sites.yaml` output.

The generator enforces an explicit acceptance criterion (`check_recovers_old`):
the new file must be a **superset** — for every `(clade, site)` in the old file
*within the `site_include_range`*, that site must still be masked for that clade
(as a whole-site mask, or now at the **mutation** level — i.e. ≥ 1 new mutation at
that position for that clade, an intentional, more precise "downgrade"). Old masks
in the UTRs are intentionally not recovered (§7) and are reported separately.

On the 2026-04-22 gisaid-and-public tree (founder-anchored), the recovery check
passes (exit 0): every old `(clade, site)` pair is preserved, with 141 whole-site →
mutation downgrades (in `22B` / `22D` / `22E`), 217 old UTR masks intentionally
dropped, and zero truly missing. `strict: true` makes any *true* shortfall a hard
error; `union_old: true` would force-include unrecovered old masks.

### Per-clade breakdown of added masks

Beyond recovering the old masks, the new scheme adds many more. The table below
counts, per Nextstrain clade, the in-range **masked positions** (a position counts
whether it is masked whole-site *or* via a mutation) in the old vs. new scheme, and
how many the new scheme **adds** (masked now, not before). This is a snapshot of
`results/recovery_summary.txt`, which the `check_recovery` rule regenerates on every
run; refresh the snapshot below from that file when the spec/tree change. Measured
on the 2026-04-22 tree:

| clade | old | new | added | | clade | old | new | added |
|---|--:|--:|--:|---|---|--:|--:|--:|
| 20B | 0 | 9 | **9** *(new clade)* | | 23B | 0 | 82 | **82** *(new clade)* |
| 20H | 9 | 9 | 0 | | 23C | 0 | 92 | **92** *(new clade)* |
| 20I | 9 | 9 | 0 | | 23D | 0 | 82 | **82** *(new clade)* |
| 20J | 9 | 9 | 0 | | 23E | 0 | 89 | **89** *(new clade)* |
| 21A | 23 | 23 | 0 | | 23F | 0 | 82 | **82** *(new clade)* |
| 21I | 23 | 23 | 0 | | 23G | 0 | 82 | **82** *(new clade)* |
| 21J | 23 | 23 | 0 | | 23H | 0 | 82 | **82** *(new clade)* |
| 21K | 63 | 64 | **1** | | 23I | 0 | 109 | **109** *(new clade)* |
| 21L | 11 | 37 | **26** | | 24A | 0 | 117 | **117** *(new clade)* |
| 21M | 0 | 20 | **20** *(new clade)* | | 24B | 0 | 118 | **118** *(new clade)* |
| 22A | 26 | 52 | **26** | | 24C | 0 | 117 | **117** *(new clade)* |
| 22B | 61 | 86 | **25** | | 24D | 0 | 117 | **117** *(new clade)* |
| 22C | 11 | 37 | **26** | | 24E | 0 | 128 | **128** *(new clade)* |
| 22D | 66 | 93 | **27** | | 24F | 0 | 117 | **117** *(new clade)* |
| 22E | 61 | 88 | **27** | | 24G | 0 | 117 | **117** *(new clade)* |
| 22F | 11 | 83 | **72** | | 24H | 0 | 117 | **117** *(new clade)* |
| 23A | 0 | 82 | **82** *(new clade)* | | 24I | 0 | 117 | **117** *(new clade)* |
| | | | | | 25A | 0 | 117 | **117** *(new clade)* |
| | | | | | 25B | 0 | 117 | **117** *(new clade)* |
| | | | | | 25C | 0 | 117 | **117** *(new clade)* |

**Totals: 406 → 2863 masked positions, +2457 added.** Of those, **2227** fall in
the **23 clades absent from the old hand-made scheme** (`20B`, `21M`, and all of
`23*`/`24*`/`25*`) and **230** are added within the **14 clades the old scheme
already had**. The four oldest clades (`20H`/`20I`/`20J`, `21A`/`21I`/`21J`) gain
nothing — the old scheme already covered them fully. Counting *distinct genomic
positions* (across all clades) rather than clade-position pairs, the scheme grows
from 161 → 273 masked positions (+112 never masked in any clade before).

The added masks come from **topological inheritance** down the tree: each
Nextstrain block accumulates the masks of *every* PANGO spec clade on its lineage
(the `# accumulated from PANGO clades: [...]` comment in each block records which),
whereas the 2022 hand-made scheme was transcribed only from the early
`maskDelta.sh` era and covered a subset of clades and sites. Two illustrative
cases:

- **`21L` (BA.2), +26.** Accumulates from `['BA.2']`. The old scheme masked only
  BA.2's recurring `11288–11296` deletion and the homoplasic `28877`/`28878` pair;
  the new scheme adds the rest of the in-range `BA.2` spec entry — two more
  deletions (`21633–21641`, `28361–28371`), three amplicon-dropout sites
  (`22786`, `22882`, `23854`), and three recurrent multi-mutation sites (`28245`,
  `28251`, `28254`). (The `BA.2` entry's UTR sites are dropped by the range filter.)
- **`22F` (XBB), +72.** Much larger because XBB sits deeper in the tree and is a
  recombinant. Its block accumulates from `['BA.2', 'XBB', 'XBB.1', 'XBB.1.5']`:
  the same 26 BA.2 whole-site masks, plus XBB's own deletion/homopolymer sites
  (`21991–21995`, `21998`; 6 positions), plus XBB's large 38-mutation `reversions`
  list (XBB has "lots of false reversions"), plus one reversion each from `XBB.1`
  (`T22317G`) and `XBB.1.5` (`T27915G`) — 26 + 6 + 40 = 72, of which the 40
  reversions are masked at the mutation level rather than whole-site.

## Updating the vendored spec

Re-download `branchSpecificMask.yml` from the UCSC kent repo
(`src/hg/utils/otto/sarscov2phylo/branchSpecificMask.yml`) at a newer commit,
update the provenance block at the top of `data/branchSpecificMask.yml`
(current: commit `2c7413b`, 2026-03-13), then regenerate as above. If the run
fails — a new clade with no founder node, an anchor cross-check violation, or an
unlisted absent representative — confirm the sample by hand in the tree and add an
approved entry to `representative_name_overrides` (renamed) or
`representative_absent_ok` (genuinely absent) in `config.yaml` — the script never
guesses. Raise `max_anchor_distance` only after confirming by hand that the founder
and representative anchors really are close.

## Caveats / accepted approximations

- **Coarsening loss (accepted).** PANGO clades finer than any Nextstrain clade
  are widened to the enclosing Nextstrain clade (over-mask). Logged per clade.
- **Reversions.** Handled at the mutation level, not blanket-masked over the site.
- **Exclusions.** Per-sample un-masking of a site cannot be expressed at
  Nextstrain granularity → skipped + logged (currently only BA.2.86 site 23599).
- **Recombinants.** The non-monophyletic `recombinant` clade has no single
  founder node and is left unmasked.
- **Founder anchoring (manual sign-off).** Masks are anchored at each clade's
  founder (MRCA); the representative (matched by exact pre-`|` name only) is a
  cross-check bounded by `max_anchor_distance`. Name overrides and absent-ok
  entries each require a hard-coded, human-reviewed config entry (see §3).
- **GISAID tree required.** Several representatives are GISAID-only; outputs
  inherit GISAID redistribution terms (the YAML itself carries no sequence data).
- **Drift.** The UCSC spec grows as new lineages are designated; the vendored copy
  is pinned to a specific commit and should be refreshed periodically.
