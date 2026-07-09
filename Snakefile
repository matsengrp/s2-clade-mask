# =============================================================================
# s2-clade-mask pipeline
#
# Staged generation of results/usher_masked_sites.yaml from the UCSC
# PANGO-lineage-keyed branchSpecificMask.yml, decomposed so each expensive stage
# is cached and independently re-runnable (the two matUtils passes especially --
# the --sample-paths pass is ~1 hour and produces a ~19 GB file).
#
#   conda activate s2-clade-mask
#   snakemake --cores 4              # full run
#   snakemake -n                     # dry run (show the DAG)
#
# All run parameters come from config.yaml, read natively below. Run from the
# repo root so the relative paths in config.yaml resolve. To skip the matUtils
# passes during development, point `sample_paths` / `clade_paths` in config.yaml
# at precomputed TSVs -- those files then satisfy the inputs and the extract
# rules never fire.
# =============================================================================

configfile: "config.yaml"

# Precomputed-paths handling: if config supplies a path, use it (the file
# satisfies the input and the corresponding extract rule is not run); otherwise
# produce it under results/.
SAMPLE_PATHS = config.get("sample_paths") or "results/sample_paths.tsv"
CLADE_PATHS  = config.get("clade_paths")  or "results/clade_paths.tsv"

NODES_REP     = config.get("nodes_representative", "results/nodes_representative.yaml")
NODES_FOUNDER = config.get("nodes_founder", "results/nodes_founder.yaml")
OUTPUT        = config["output"]
OUTPUT_PANGO  = config.get("output_pango", "results/usher_masked_sites_pango.yaml")
COMPARISON    = config.get("comparison_report", "results/anchoring_comparison.txt")
DISTANCES     = config.get("anchor_distances", "results/anchor_distances.yaml")
ANCHOR_CHECK  = config.get("anchor_check_report", "results/anchor_check.txt")
RECOVERY      = config.get("recovery_report", "results/recovery_report.txt")
RECOVERY_SUMMARY = config.get("recovery_summary", "results/recovery_summary.txt")
BRANCH_SPEC   = config["branch_spec"]
METADATA      = config["metadata"]
OLD_YAML      = config.get("old_yaml")

CONFIG = workflow.configfiles[0] if workflow.configfiles else "config.yaml"


rule all:
    input:
        OUTPUT,
        OUTPUT_PANGO,
        COMPARISON,
        ANCHOR_CHECK,
        RECOVERY,
        RECOVERY_SUMMARY,


rule extract_sample_paths:
    """matUtils extract --sample-paths (root->tip node chain per sample; ~19 GB)."""
    input:
        tree=config["tree"] if config.get("tree") else [],
    output:
        SAMPLE_PATHS,
    shell:
        "python scripts/extract_paths.py --config {CONFIG} --kind sample --out {output}"


rule extract_clade_paths:
    """matUtils extract --clade-paths ({clade: founder node id})."""
    input:
        tree=config["tree"] if config.get("tree") else [],
    output:
        CLADE_PATHS,
    shell:
        "python scripts/extract_paths.py --config {CONFIG} --kind clade --out {output}"


rule resolve_representative:
    """Mode 1 (string match): PANGO clade -> representative sample node."""
    input:
        sample_paths=SAMPLE_PATHS,
        branch_spec=BRANCH_SPEC,
    output:
        NODES_REP,
    shell:
        "python scripts/resolve_nodes.py --config {CONFIG} --mode representative "
        "--sample-paths {input.sample_paths} --out {output}"


rule resolve_founder:
    """Mode 2 (founder): PANGO clade -> annotated founder (MRCA) node."""
    input:
        clade_paths=CLADE_PATHS,
        branch_spec=BRANCH_SPEC,
    output:
        NODES_FOUNDER,
    shell:
        "python scripts/resolve_nodes.py --config {CONFIG} --mode founder "
        "--clade-paths {input.clade_paths} --out {output}"


rule compare_and_build:
    """Emit the founder-anchored Nextstrain + PANGO YAMLs and the anchor table."""
    input:
        nodes_representative=NODES_REP,
        nodes_founder=NODES_FOUNDER,
        sample_paths=SAMPLE_PATHS,
        clade_paths=CLADE_PATHS,
        metadata=METADATA,
        branch_spec=BRANCH_SPEC,
    output:
        yaml=OUTPUT,
        pango=OUTPUT_PANGO,
        comparison=COMPARISON,
        distances=DISTANCES,
    shell:
        "python scripts/compare_and_build.py --config {CONFIG} "
        "--nodes-representative {input.nodes_representative} "
        "--nodes-founder {input.nodes_founder} "
        "--sample-paths {input.sample_paths} --clade-paths {input.clade_paths} "
        "--out-yaml {output.yaml} --out-pango {output.pango} "
        "--out-comparison {output.comparison} --out-distances {output.distances}"


rule check_anchors:
    """Fail if any founder anchor is too far (mutations) from its representative."""
    input:
        distances=DISTANCES,
    output:
        ANCHOR_CHECK,
    shell:
        "python scripts/check_anchors.py --config {CONFIG} "
        "--distances {input.distances} --out {output}"


rule check_recovery:
    """Verify the Nextstrain-clade YAML is a superset of the vendored old scheme."""
    input:
        yaml=OUTPUT,
        old=OLD_YAML if OLD_YAML else [],
    output:
        report=RECOVERY,
        summary=RECOVERY_SUMMARY,
    shell:
        "python scripts/check_recovery.py --config {CONFIG} "
        "--yaml {input.yaml} --out {output.report} --out-summary {output.summary}"
