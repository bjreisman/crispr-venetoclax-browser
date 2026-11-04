# Methods and interpretation

This is an exploratory browser of heterogeneous published results, not a validated pooled effect-size meta-analysis. The initial public version preserves the existing score transformations and analysis results.

## Included evidence

The core cache contains ten studies. It combines survival-selection knockout screens, ORF overexpression and CRISPR activation screens, BH3-priming sorting, and baseline essentiality. Treatment arms include venetoclax, other BH3 mimetics, and combinations. Compare source-native results within their experimental context.

`data/metadata/included_screen_inventory.csv` records biological screen units and source links. `experiment_to_screen_mapping.csv` maps displayed experiments to those units. Replicates, timepoints, selection tests, and derived summaries do not establish independent biological replication. The app's current recurrence table still counts dataset presence and experiment labels; it does not use these registries to calculate replicated directional support.

## Current limitations

- The significant/annotated flag can use nominal p <= 0.05 when FDR is available but exceeds 0.05. A positive flag therefore does not imply adjusted significance. Missing FDR values pass the current FDR filter.
- MAGeCK positive- and negative-selection arms both report the signed gene-level LFC; the arm determines only which p-value, FDR, and rank are shown. Releases before v0.1.1 imposed the arm's sign on the effect magnitude, so enriched genes appeared depleted in negative-selection arms.
- Raw scores use different units. Neither raw-effect heatmap colors nor a common magnitude threshold establish comparable biological effects across studies. Centering an experiment to a z-score can change an effect's sign.
- Eide2023 guide-count summaries scale each sample to the median library size and center each experiment on its median gene (the library has no non-targeting guides). They are still exploratory treatment/control ratios without replicate-aware scoring, guide/control QC, p-values, or FDRs.
- Patient expression/response correlations may include repeated specimens and measurements from the same subject; reported row counts are not independent-subject counts.
- Gene matching normalizes symbols and uses curated aliases; it does not perform formal human/mouse orthology mapping. TP53 and mouse Trp53 remain separate keys. The four mouse A1 paralogs (Bcl2a1a/b/c/d) map to human BCL2A1; in the Gene Explorer screen view, merged paralogs keep their most extreme value.
- Paper Explorer previews at most 10,000 rows. Its table CSV can be a subset; use the separate Download View control for the full filtered view.

These limitations affect interpretation and target prioritization. They do not prevent sharing the software as a work in progress. The next scientific milestones are in [PUBLICATION_PLAN.md](PUBLICATION_PLAN.md).

## Gene Explorer screen summary

The default Gene Explorer view shows one representative result per independent screen. `R/screen_views.R` defines the representative experiments, using the parent screens in `data/metadata/experiment_to_screen_mapping.csv`:

- Repeat timepoints use the latest timepoint. Guièze2019 STARS uses replicate 1.
- MAGeCK screens list both selection arms. Each gene takes the arm that matches its LFC sign, so the displayed FDR is the directional one.
- Derived aggregates, hit subsets, and Lin2019 essentiality scores are excluded.
- Screens are grouped as venetoclax in AML, venetoclax in lymphoid models, and venetoclax combinations (shown by default), plus other BH3 mimetics and CRISPRa, ORF, and BH3-priming assays (optional).

The axis is a signed within-screen percentile, from −1 (most depleted) to +1 (most enriched), computed among all genes in that screen. Lin2019 percentiles are relative to its metabolic library. Guièze2019 STARS, Portelinha2024, and the Deng2022 CRISPRa tables list only selected genes, so listed genes are shown without a percentile.

The per-gene summary counts venetoclax screens in which the gene falls in the top or bottom 10%, and how many of those also have FDR ≤ 0.05. Hit-list screens count as tested for every gene: a listed gene with FDR ≤ 0.05 counts as resistance, and an unlisted gene counts as no signal. This is a descriptive summary, not a meta-analysis. Percentiles are not comparable between ORF or BH3-priming assays and knockout survival screens.

## Study-specific processing

Guièze2019 includes STARS knockout results and ORF overexpression. The ORF browser uses the strongest venetoclax-enriched construct per gene; it is a different perturbation from knockout.

He2026 includes Figure 1 combined MAGeCK results with p-values and P1-R/P2-R sgRNA RPM-derived LFC summaries without p-values/FDRs. The LFC is CytoC-high/CytoC-low; the current browser interprets negative values as CytoC-low enrichment.

Nakao2023 uses DrugZ normalized Z as its drug-vs-DMSO effect, with directional p-values/FDRs. MAGeCK MLE results remain audit fields. DrugZ v1.1.0.2 build 116 uses 10-million-read normalization, pseudocount 5, and half-window 500. MAGeCK 0.5.9.5 uses median normalization. Non-targeting and miRNA pseudo-genes remain in the analysis but are omitted from browser rows. The original local analysis manifest records source/output SHA-256 hashes and commands; full source workbooks and intermediate results are not needed to run the prepared browser.

Deng2022's genome-wide CRISPR activation screens use the published Supplementary Data 2 tables: significantly enriched sgRNAs from edgeR/voom, IC80 venetoclax (30 nM for line #214, 15 nM for #216) vs DMSO. They are imported as hit lists, with each gene key keeping its strongest guide (lowest FDR, then largest logFC). The raw counts in Supplementary Data 3 are heavily bottlenecked: fewer than 40% of guides are detected even before treatment, and about 1% after. They are not reanalyzed into genome-wide scores. CRISPRa is a gain-of-function perturbation, so enrichment means that activating the gene confers resistance; it is grouped with the ORF screen, not with the knockouts. The paper's lymphoma vs B-ALL RNA-seq comparison is not imported: it compares tumour types (8 vs 2 samples, sequenced on different platforms) rather than measuring venetoclax response, and it served to validate the lymphoma model.

DepMap Public 26Q1 contributes baseline CRISPR gene effects for MOLM-13, THP-1, and OCI-LY19. OCI-Ly1 has no direct match in this layer. Baseline dependency is context, not a treatment-response effect or an established therapeutic window.

Aliases group ADSS with ADSS2, ADSS1 with ADSSL1, and TOMM70A with TOMM70. ADSSL1 remains distinct from ADSS2.

## Additional-screen catalog

The 31 October 2026 audit records 29 candidate entries in `data/metadata/additional_screen_candidates.csv`, with scope and limitations in `additional_screen_search_notes.txt`. These are candidates and supporting references, not additional imported scores. EndNote record numbers are local evidence locators; public source links are supplied in the catalog.
