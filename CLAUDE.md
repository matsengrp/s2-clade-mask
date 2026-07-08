# CLAUDE.md

Guidance for Claude Code (claude.ai/code) and other contributors working in this
repository.

## Project overview

`s2-clade-mask` is a **standalone, single-purpose** repository. Its only job is to
generate `results/usher_masked_sites.yaml` — the per-**Nextstrain-clade**
site/mutation masking file consumed by the
[`mat-pcp`](https://github.com/matsengrp/mat-pcp) SARS-CoV-2 pipeline — from the
UCSC PANGO-lineage-keyed branch-specific masking spec
(`data/branchSpecificMask.yml`).

When the masks are refreshed here, the resulting YAML is **copied into `mat-pcp`**
at `data/usher_masked_sites.yaml`. That copy is the only thing `mat-pcp` needs
from this repo; generation (which needs the full ~8M-sample MAT + `matUtils`) is
fully decoupled from application.

Read `README.md` for the full algorithm write-up; this file covers how to work in
the repo without re-deriving that.

## Layout

```
config.yaml                            # all run parameters (paths + options)
environment.yml                        # conda env: usher (matUtils), snakemake, python, pyyaml, gh
Snakefile                              # staged pipeline (primary run path)
scripts/s2clademask.py                 # shared library: ALL generation logic lives here
scripts/extract_paths.py               # rule: one matUtils extract pass (--sample-paths/--clade-paths)
scripts/resolve_nodes.py               # rule: resolve clades in one anchoring mode
scripts/compare_and_build.py           # rule: compare both modes (report) + emit final YAML
scripts/check_recovery.py              # rule: old-scheme superset recovery check
scripts/generate_usher_masked_sites.py # one-shot wrapper (same lib, single process)
data/branchSpecificMask.yml            # vendored UCSC input spec (pinned commit)
data/usher_masked_sites_old_scheme.yml # vendored jbloomlab old-scheme baseline (pinned commit)
results/usher_masked_sites.yaml        # generated output (tracked)
```

The generation logic lives in **one place** — `scripts/s2clademask.py`. Both the
Snakemake rule scripts and the one-shot wrapper import it; **never** duplicate
logic into a rule script. The pipeline and the wrapper must produce an identical
`results/usher_masked_sites.yaml` because they call the same functions.

## Running

```bash
conda env create -f environment.yml
conda activate s2-clade-mask
# from the repo root, so relative paths in config.yaml resolve:
snakemake --cores 4                                        # primary path (staged, cached)
# or the equivalent one-shot process:
python scripts/generate_usher_masked_sites.py --config config.yaml
```

The Snakemake pipeline is the primary path: `extract_sample_paths` /
`extract_clade_paths` (the two matUtils passes) → `resolve_representative` (Mode 1,
string match) + `resolve_founder` (Mode 2, MRCA) → `compare_and_build` (writes the
report-only `anchoring_comparison.txt` and the final YAML) → `check_recovery`.
Each stage is cached, so troubleshooting doesn't re-pay for the ~1 h / ~19 GB
`--sample-paths` extract. During development, set `sample_paths` / `clade_paths`
in `config.yaml` to precomputed TSVs so the extract rules are skipped.

All run parameters come from `config.yaml`. Both the pipeline and the wrapper run
against the **full gisaid-and-public** MAT — never the public-only tree (missing
several GISAID-only representatives) or a filtered per-run tree.

## Design invariants — do not break these

1. **Representative-primary resolution.** Each PANGO clade is anchored at its
   representative sample's node (matching upstream `matUtils mask`), not the clade
   MRCA. Founder nodes are only an approved per-clade fallback.
2. **Exact matching only.** Representatives match the tree sample name by exact
   equality on the field **before the first `|`**. Never add fuzzy/substring
   matching. An unresolved representative with no approved entry must remain a
   **hard error**, never a silent guess.
3. **Every deviation is human-approved.** The `representative_name_overrides`
   (renamed samples) and `representative_fallbacks` (absent samples) tables live in
   `config.yaml` and are validated at load time. The config is version-controlled,
   so entries are reviewed exactly like code. Add an entry — with a comment naming
   the clade and why — only after confirming the sample by hand; never add fuzzy
   matching to avoid one.
4. **Superset recovery.** The output must remain a superset of the old hand-made
   scheme within `site_include_range` (checked via `old_yaml` + `strict`). `old_yaml`
   always points at the vendored `data/usher_masked_sites_old_scheme.yml` — the
   original jbloomlab/SARS2-mut-fitness file, pinned at a fixed commit — never at a
   previous run's `results/usher_masked_sites.yaml`. A true shortfall is a bug to
   investigate, not to accept.
5. **UTR masks dropped.** Masks outside `[site_include_start, site_include_end]`
   are omitted (mat-pcp discards those positions at runtime anyway).

## Editing conventions

- **All generation logic → `scripts/s2clademask.py`.** The rule scripts
  (`extract_paths.py`, `resolve_nodes.py`, `compare_and_build.py`,
  `check_recovery.py`) and the one-shot `generate_usher_masked_sites.py` are thin
  entrypoints that import the library; do not reimplement logic in them. Adding a
  Snakemake rule means adding a thin script + a rule, not a second copy of an
  algorithm.
- **All run inputs** — tree/metadata paths, output, ranges, strict, the pipeline
  intermediate/report paths, *and* the reviewed masking decisions
  (`representative_name_overrides` / `representative_fallbacks`) → `config.yaml`.
  Keep the schema in `CONFIG_DEFAULTS` / `REQUIRED_CONFIG` (in `s2clademask.py`)
  and the `validate_resolution_overrides` checks in sync with `config.yaml` and
  the README key table.
- **Never hand-edit** `results/usher_masked_sites.yaml` (regenerate it),
  `data/branchSpecificMask.yml`, or `data/usher_masked_sites_old_scheme.yml`
  (re-vendor them — see below).
- Keep `README.md` and this file **in sync** with any change to layout, the config
  schema, the resolution logic, the pipeline/rule structure, or the run procedure.

## Updating the vendored spec

Re-download `branchSpecificMask.yml` from the UCSC kent repo
(`src/hg/utils/otto/sarscov2phylo/branchSpecificMask.yml`) at a newer commit,
update the provenance block at the top of `data/branchSpecificMask.yml`, then
regenerate. If a representative was renamed or dropped upstream, the run fails
loudly — confirm the sample by hand in the tree and add an approved
`representative_name_overrides` (renamed) or `representative_fallbacks` (absent)
entry to `config.yaml`. After regenerating, re-measure the count tables in the README.

`data/usher_masked_sites_old_scheme.yml` (the `old_yaml` baseline) is a fixed
historical reference, not something to refresh routinely — it should only ever
change if the acceptance criterion itself is revisited (e.g. a different upstream
commit is agreed to be the authoritative baseline). If that happens, re-download
from the URL in its provenance header and update the block the same way.

## Sanity checks

- `python scripts/generate_usher_masked_sites.py --help` (imports resolve).
- `load_config` rejects unknown keys and missing required keys.
- A full run with `strict: true` exits 0 only if old-scheme recovery is clean and
  every representative resolves (or has an approved override/fallback).
