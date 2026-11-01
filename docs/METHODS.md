# Methods and interpretation

This is an exploratory browser of heterogeneous published results, not a validated pooled effect-size meta-analysis. The initial public version preserves the existing score transformations and analysis results.

## Included evidence

The core cache contains ten studies. It combines survival-selection knockout screens, ORF overexpression, BH3-priming sorting, baseline essentiality, and contextual RNA-seq differential expression. Treatment arms include venetoclax, other BH3 mimetics, and combinations. Compare source-native results within their experimental context.

`data/metadata/included_screen_inventory.csv` records biological screen units and source links. `experiment_to_screen_mapping.csv` maps displayed experiments to those units. Replicates, timepoints, selection tests, and derived summaries do not establish independent biological replication. The app's current recurrence table still counts dataset presence and experiment labels; it does not use these registries to calculate replicated directional support.

## Current limitations

- The significant/annotated flag can use nominal p <= 0.05 when FDR is available but exceeds 0.05. A positive flag therefore does not imply adjusted significance. Missing FDR values pass the current FDR filter.
- Some MAGeCK importers impose the positive/negative selection-arm sign on effect magnitude. Inspect source contrasts before interpreting biological direction.
- Raw scores use different units. Neither raw-effect heatmap colors nor a common magnitude threshold establish comparable biological effects across studies. Centering an experiment to a z-score can change an effect's sign.
- Eide2023 guide-count summaries use exploratory treatment/control ratios without a completed replicate-aware normalization and QC analysis.
- Patient expression/response correlations may include repeated specimens and measurements from the same subject; reported row counts are not independent-subject counts.
- Gene matching normalizes symbols and uses curated aliases; it does not perform formal human/mouse orthology mapping. TP53 and mouse Trp53 remain separate keys.
- Paper Explorer previews at most 10,000 rows. Its table CSV can be a subset; use the separate Download View control for the full filtered view.

These limitations affect interpretation and target prioritization. They do not prevent sharing the software as a work in progress. The next scientific milestones are in [PUBLICATION_PLAN.md](PUBLICATION_PLAN.md).

## Study-specific processing

Guièze2019 includes STARS knockout results and ORF overexpression. The ORF browser uses the strongest venetoclax-enriched construct per gene; it is a different perturbation from knockout.

He2026 includes Figure 1 combined MAGeCK results with p-values and P1-R/P2-R sgRNA RPM-derived LFC summaries without p-values/FDRs. The LFC is CytoC-high/CytoC-low; the current browser interprets negative values as CytoC-low enrichment.

Nakao2023 uses DrugZ normalized Z as its drug-vs-DMSO effect, with directional p-values/FDRs. MAGeCK MLE results remain audit fields. DrugZ v1.1.0.2 build 116 uses 10-million-read normalization, pseudocount 5, and half-window 500. MAGeCK 0.5.9.5 uses median normalization. Non-targeting and miRNA pseudo-genes remain in the analysis but are omitted from browser rows. The original local analysis manifest records source/output SHA-256 hashes and commands; full source workbooks and intermediate results are not needed to run the prepared browser.

Deng2022 contributes contextual CRISPRa-associated RNA-seq differential expression, rather than directly equivalent knockout scores.

DepMap Public 26Q1 contributes baseline CRISPR gene effects for MOLM-13, THP-1, and OCI-LY19. OCI-Ly1 has no direct match in this layer. Baseline dependency is context, not a treatment-response effect or an established therapeutic window.

Aliases group ADSS with ADSS2, ADSS1 with ADSSL1, and TOMM70A with TOMM70. ADSSL1 remains distinct from ADSS2.

## Additional-screen catalog

The 31 October 2026 audit records 29 candidate entries in `data/metadata/additional_screen_candidates.csv`, with scope and limitations in `additional_screen_search_notes.txt`. These are candidates and supporting references, not additional imported scores. EndNote record numbers are local evidence locators; public source links are supplied in the catalog.
