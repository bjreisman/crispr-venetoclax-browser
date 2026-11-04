# Methods and interpretation

This is an exploratory browser of heterogeneous published results, not a validated pooled effect-size meta-analysis. The initial public version preserves the existing score transformations and analysis results.

## Included evidence

The core cache contains sixteen studies. It combines survival-selection knockout screens, ORF overexpression and CRISPR activation screens, BH3-priming sorting, and baseline essentiality. Treatment arms include venetoclax, other BH3 mimetics, and combinations. Compare source-native results within their experimental context.

`data/metadata/included_screen_inventory.csv` records biological screen units and source links. `experiment_to_screen_mapping.csv` maps displayed experiments to those units. Replicates, timepoints, selection tests, and derived summaries do not establish independent biological replication. The app's current recurrence table still counts dataset presence and experiment labels; it does not use these registries to calculate replicated directional support.

## Current limitations

- The significant/annotated flag can use nominal p <= 0.05 when FDR is available but exceeds 0.05. A positive flag therefore does not imply adjusted significance. Missing FDR values pass the current FDR filter.
- MAGeCK positive- and negative-selection arms both report the signed gene-level LFC; the arm determines only which p-value, FDR, and rank are shown. Releases before v0.1.2 imposed the arm's sign on the effect magnitude, so enriched genes appeared depleted in negative-selection arms.
- Raw scores use different units. Neither raw-effect heatmap colors nor a common magnitude threshold establish comparable biological effects across studies. Centering an experiment to a z-score can change an effect's sign.
- Eide2023 guide-count summaries scale each sample to the median library size and center each experiment on its median gene (the library has no non-targeting guides). They are still exploratory treatment/control ratios without replicate-aware scoring, guide/control QC, p-values, or FDRs.
- Patient expression/response correlations may include repeated specimens and measurements from the same subject; reported row counts are not independent-subject counts.
- Gene matching normalizes symbols and uses curated aliases; it does not perform formal human/mouse orthology mapping. TP53 and mouse Trp53 remain separate keys. The four mouse A1 paralogs (Bcl2a1a/b/c/d) map to human BCL2A1; in the Gene Explorer screen view, merged paralogs keep their most extreme value.
- The Paper Explorer full table previews at most 10,000 rows. Its table CSV can be a subset; use the separate Download View control for every row of the selected result.

These limitations affect interpretation and target prioritization. They do not prevent sharing the software as a work in progress. The next scientific milestones are in [PUBLICATION_PLAN.md](PUBLICATION_PLAN.md).

## Paper Explorer

The Paper Explorer shows one result at a time: one comparison, with MAGeCK positive- and negative-selection arms merged as in the Gene Explorer. Results are grouped by the screen group of their parent screen; ★ marks the representative result used in the Gene Explorer, and Lin2019 essentiality scores and the Stewart2026 palbociclib-only arm sit under "Other tables".

- The plot orders every gene by within-screen percentile (x) and draws it on the chosen metric (y): robust z (default), signed −log₁₀ FDR as published, or the raw effect as published. Hit tables without genome-wide scores order genes by their published effect.
- The significance counts use FDR ≤ 0.05 or the paper's own threshold. Stewart2026 venetoclax and combination arms show "—" because their published FDRs come from unnormalized values.
- "Other venetoclax screens" counts the venetoclax screens from other papers (the Gene Explorer's default groups) where the gene is in the same tail: top or bottom 10% of the screen, or listed in a published resistance hit table. Screens from the same paper are excluded so that timepoints and replicates do not count as replication.

## Gene Explorer screen summary

The default Gene Explorer view shows one representative result per independent screen. `R/screen_views.R` defines the representative experiments, using the parent screens in `data/metadata/experiment_to_screen_mapping.csv`:

- Repeat timepoints use the latest timepoint (Sharon2019 uses day 29). Guièze2019 STARS uses replicate 1.
- MAGeCK screens list both selection arms. Each gene takes the arm that matches its LFC sign, so the displayed FDR is the directional one.
- Derived aggregates, hit subsets, and Lin2019 essentiality scores are excluded.
- Screens are grouped as venetoclax in AML, venetoclax in lymphoid models, and venetoclax combinations (shown by default), plus other BH3 mimetics and CRISPRa, ORF, and BH3-priming assays (optional).

The axis is a signed within-screen percentile, from −1 (most depleted) to +1 (most enriched), computed among all genes in that screen. Lin2019 percentiles are relative to its metabolic library. Guièze2019 STARS, Portelinha2024, and the Deng2022 and Diepstraten2026 CRISPRa tables list only selected genes, so listed genes are shown without a percentile.

Screens with more than one result (repeat timepoints, replicates, doses, or derived tables such as Lin2019's average and Nechiporuk2019's top-hits table) are marked ▸. Clicking one of the paper's dots, or choosing it under Expand papers, shows each result as an indented row, scored within its own experiment. Results are grouped by their parent screen in `data/metadata/experiment_to_screen_mapping.csv`; MAGeCK selection arms of one comparison form one row, and baseline essentiality scores are not shown. Expanded rows are for inspection only and do not enter the per-gene summary.

The per-gene summary counts venetoclax screens in which the gene falls in the top or bottom 10%, and how many of those also have FDR ≤ 0.05. Hit-list screens count as tested for every gene: a listed gene with FDR ≤ 0.05 counts as resistance, and an unlisted gene counts as no signal. This is a descriptive summary, not a meta-analysis. Percentiles are not comparable between ORF or BH3-priming assays and knockout survival screens.

## Study-specific processing

Guièze2019 includes STARS knockout results and ORF overexpression. The ORF browser uses the strongest venetoclax-enriched construct per gene; it is a different perturbation from knockout.

He2026 includes Figure 1 combined MAGeCK results with p-values and P1-R/P2-R sgRNA RPM-derived LFC summaries without p-values/FDRs. The LFC is CytoC-high/CytoC-low; the current browser interprets negative values as CytoC-low enrichment.

Nakao2023 uses DrugZ normalized Z as its drug-vs-DMSO effect, with directional p-values/FDRs. MAGeCK MLE results remain audit fields. DrugZ v1.1.0.2 build 116 uses 10-million-read normalization, pseudocount 5, and half-window 500. MAGeCK 0.5.9.5 uses median normalization. Non-targeting and miRNA pseudo-genes remain in the analysis but are omitted from browser rows. The original local analysis manifest records source/output SHA-256 hashes and commands; full source workbooks and intermediate results are not needed to run the prepared browser.

Wang2023 contributes two MOLM-13 venetoclax knockout screens from Table S1: the Brunello genome-wide library (25 nM, day 20) and a 492-gene RNA-binding-protein domain library. Values are the published CRISPR scores (mean sgRNA log2 fold change, venetoclax vs DMSO); no p-values or FDRs were published. The other drugs in the table are not imported. The dose and timepoint of the RNA-binding-protein screen are not stated in the methods.

Sharon2019 is a genome-wide TKO v1 knockout screen in MOLM13-R1, a MOLM-13 clone with acquired venetoclax resistance (400 nM, days 8, 16, 22 and 29), from the authors' MAGeCK workbook (BioGRID ORCS Dataset 406). Both selection arms keep the signed LFC; control targets (EGFP, LacZ, luciferase, chr10 promoter/random) are dropped. The screen view uses day 29; the paper's Figure 1B ranks day 8. Its genome-wide correlation with the naive MOLM-13 screens is near zero, and few genes pass FDR ≤ 0.05.

Yokoyama2026 contributes genome-wide knockout screens in two peripheral T-cell lymphoma lines, ED40515(−) and DEL, each under venetoclax, navitoclax and the MCL1 inhibitor ABBV-467 vs DMSO (Supplementary Table 1, MAGeCK LFC and p-values; no FDR). Venetoclax arms are lymphoid venetoclax screens; navitoclax and ABBV-467 are other BH3 mimetics. Twenty-four Excel-date-corrupted gene IDs were restored from the RefSeq accessions in the sgRNA table; two ambiguous IDs (1-Mar, 2-Mar) were dropped. The ED40515(−) venetoclax arm shows weak positive controls (BAX LFC ≈ 0.06).

Stewart2026 contributes OCI-AML2 Yusa-library knockout screens (day 21 vs DMSO) for venetoclax 200 nM, venetoclax + palbociclib 1 µM, and palbociclib alone, using the authors' edgeR gene summaries (Table S4 `mid_lfc`, `mid_pval`, `mid_fdr`). The palbociclib arm is imported but not shown in the screen view. The authors state that the venetoclax and combination contrasts could not be normalized: a few resistant clones dominated the libraries, so the published values are shifted to a median near −6 and their FDRs flag 83–98% of genes. For these two arms the browser median-centers the effect and does not use the published p-values or FDRs for significance calls, which would otherwise add a significant call to about 16,000 genes in the cross-study views. The published values are retained (`published_effect_value`, `published_p_value`, `published_fdr`) and shown by the Gene Explorer's "Raw effect (as published)" and "Signed −log₁₀ FDR" metrics and in hover text. Rank-based metrics (percentile, robust z-score) are unaffected by centering; sensitizers are compressed near the floor of the distribution, so depletion calls are weak in these screens. GEO GSE291338 counts confirm these screens are independent of Eide2023.

BrimEdwards2026 contributes an OCI-AML2 Yusa-library knockout screen with the dual BCL2/BCL-XL inhibitor AZD4320 (30 nM vs DMSO, day 14), grouped with other BH3 mimetics. The paper's gene-level table (Supplemental Data 1) is not in the published supplement, so `scripts/build_brimedwards2026.R` re-derives gene results from the GEO GSE294240 counts with the authors' edgeR pipeline (sgRNA filters, TMM, exact test, median-p sgRNA per gene). The pipeline reproduces Stewart2026's published palbociclib sheets and recovers all Tier 1 genes labelled in Figure 2A. The rebuild needs edgeR (Bioconductor), which is not in `renv.lock`; the browser itself only reads the derived table. GSE294240 has no venetoclax arm, and its samples are not shared with Eide2023 or Stewart2026.

Diepstraten2026 contributes Partita CRISPR activation screens (10 sgRNAs per TSS, five sublibraries) in the same DHL214 and DHL216 lymphoma lineage as Deng2022, venetoclax IC50 and IC80 vs DMSO. The authors ran no statistical test for these screens, and the GEO GSE296642 counts are extremely bottlenecked (0.3–0.4% of sgRNAs detected after venetoclax; replicate Spearman about 0.07). They are imported as hit lists reanalysed from the counts: an sgRNA is a hit if it reaches the authors' Figure 3B threshold of ≥ 5% read frequency in a venetoclax replicate and is higher than under DMSO; each gene keeps its strongest sgRNA. Many hits rest on a single replicate (`n_replicates_ge5pct`). p-values and FDRs are left blank.

Deng2022's genome-wide CRISPR activation screens use the published Supplementary Data 2 tables: significantly enriched sgRNAs from edgeR/voom, IC80 venetoclax (30 nM for line #214, 15 nM for #216) vs DMSO. They are imported as hit lists, with each gene key keeping its strongest guide (lowest FDR, then largest logFC). The raw counts in Supplementary Data 3 are heavily bottlenecked: fewer than 40% of guides are detected even before treatment, and about 1% after. They are not reanalyzed into genome-wide scores. CRISPRa is a gain-of-function perturbation, so enrichment means that activating the gene confers resistance; it is grouped with the ORF screen, not with the knockouts. The paper's lymphoma vs B-ALL RNA-seq comparison is not imported: it compares tumour types (8 vs 2 samples, sequenced on different platforms) rather than measuring venetoclax response, and it served to validate the lymphoma model.

DepMap Public 26Q1 contributes baseline CRISPR gene effects for MOLM-13, THP-1, and OCI-LY19. OCI-Ly1 has no direct match in this layer. Baseline dependency is context, not a treatment-response effect or an established therapeutic window.

Aliases group ADSS with ADSS2, ADSS1 with ADSSL1, and TOMM70A with TOMM70. ADSSL1 remains distinct from ADSS2.

## Additional-screen catalog

The 31 October 2026 audit records 29 candidate entries in `data/metadata/additional_screen_candidates.csv`, with scope and limitations in `additional_screen_search_notes.txt`. These are candidates and supporting references, not additional imported scores. EndNote record numbers are local evidence locators; public source links are supplied in the catalog.
