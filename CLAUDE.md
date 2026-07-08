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
environment.yml                        # conda env: usher (matUtils), python, pyyaml
scripts/generate_usher_masked_sites.py # the generator
data/branchSpecificMask.yml            # vendored UCSC input spec (pinned commit)
results/usher_masked_sites.yaml        # generated output (tracked)
```

## Running

```bash
conda env create -f environment.yml
conda activate s2-clade-mask
# from the repo root, so relative paths in config.yaml resolve:
python scripts/generate_usher_masked_sites.py --config config.yaml
```

All run parameters come from `config.yaml` (a single `--config` flag is the only
CLI argument). Run against the **full gisaid-and-public** MAT — never the
public-only tree (missing several GISAID-only representatives) or a filtered
per-run tree.

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
   scheme within `site_include_range` (checked via `old_yaml` + `strict`). A true
   shortfall is a bug to investigate, not to accept.
5. **UTR masks dropped.** Masks outside `[site_include_start, site_include_end]`
   are omitted (mat-pcp discards those positions at runtime anyway).

## Editing conventions

- **All run inputs** — tree/metadata paths, output, ranges, strict, *and* the
  reviewed masking decisions (`representative_name_overrides` /
  `representative_fallbacks`) → `config.yaml`. Keep the schema in `CONFIG_DEFAULTS`
  / `REQUIRED_CONFIG` and the `validate_resolution_overrides` checks in sync with
  `config.yaml` and the README key table.
- **Never hand-edit** `results/usher_masked_sites.yaml` (regenerate it) or
  `data/branchSpecificMask.yml` (re-vendor it — see below).
- Keep `README.md` and this file **in sync** with any change to layout, the config
  schema, the resolution logic, or the run procedure.

## Updating the vendored spec

Re-download `branchSpecificMask.yml` from the UCSC kent repo
(`src/hg/utils/otto/sarscov2phylo/branchSpecificMask.yml`) at a newer commit,
update the provenance block at the top of `data/branchSpecificMask.yml`, then
regenerate. If a representative was renamed or dropped upstream, the run fails
loudly — confirm the sample by hand in the tree and add an approved
`representative_name_overrides` (renamed) or `representative_fallbacks` (absent)
entry to `config.yaml`. After regenerating, re-measure the count tables in the README.

## Sanity checks

- `python scripts/generate_usher_masked_sites.py --help` (imports resolve).
- `load_config` rejects unknown keys and missing required keys.
- A full run with `strict: true` exits 0 only if old-scheme recovery is clean and
  every representative resolves (or has an approved override/fallback).
