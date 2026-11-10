# Data sources and attribution

Original project code is MIT licensed. The numerical data, source publications, and third-party tools retain their own upstream terms. Do not interpret the software license as a new license for these materials.

## Public source datasets

The prepared caches are derived from public study supplements and repositories. The [source manifest](../data/metadata/source_manifest.csv) records required local paths, public source links, and input roles. The [screen inventory](../data/metadata/included_screen_inventory.csv) and [experiment mapping](../data/metadata/experiment_to_screen_mapping.csv) provide more specific evidence locators.

- [Chen2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6606342/)
- [Deng2022](https://doi.org/10.1038/s41467-022-32485-9): CRISPRa enriched-sgRNA tables from Supplementary Data 2 ([PMC9374748](https://pmc.ncbi.nlm.nih.gov/articles/PMC9374748/))
- [BrimEdwards2026](https://doi.org/10.3324/haematol.2025.287972): counts in [GSE294240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294240), re-analysed with the authors' edgeR pipeline
- [Diepstraten2026](https://doi.org/10.1126/sciadv.aec0722): counts in [GSE296642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296642), re-analysed as hit lists
- [Eide2023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216087)
- [Glytsou2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10330144/)
- [Guièze2019](https://jourdainlab.org/wp-content/uploads/2023/03/Guieze_CancerCell_2019.pdf)
- [He2026](https://www.nature.com/articles/s43018-026-01184-5)
- [Lin2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6506362/)
- [Nakao2023](https://www.nature.com/articles/s41375-023-01879-z)
- [Nechiporuk2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6606338/)
- [Portelinha2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11380151/)
- [Sharon2019](https://doi.org/10.1126/scitranslmed.aax2863): MAGeCK workbook from [BioGRID ORCS Dataset 406](https://orcs.thebiogrid.org/Dataset/406)
- [Stewart2026](https://doi.org/10.1016/j.xcrm.2025.102526): Table S4 ([PMC12866115](https://pmc.ncbi.nlm.nih.gov/articles/PMC12866115/)); counts in [GSE291338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291338)
- [Wang2023](https://doi.org/10.1016/j.ccell.2022.12.002): Table S1 ([PMC9839614](https://pmc.ncbi.nlm.nih.gov/articles/PMC9839614/))
- [Yokoyama2026](https://doi.org/10.1038/s41388-026-03970-y): Supplementary Tables 1 and 2

## External context

- **BeatAML2:** [Bottomly and colleagues’ data repository](https://github.com/biodev/beataml2.0_data). The authors’ [download page](https://biodev.github.io/BeatAML2/) explicitly licenses these processed data under [CC BY 4.0](https://github.com/biodev/beataml2.0_data/blob/main/LICENSE.txt). The cache selects venetoclax dose-response rows, computes viability AUCs, aggregates expression by normalized gene symbol, and calculates curated gene-set scores. These are derived/modified data; cite the original BeatAML2 manuscript and data repository. Public dbGaP-labelled subject/sample identifiers are retained to connect response and expression measurements.
- **DepMap Public 26Q1, Broad Institute:** [DepMap downloads](https://depmap.org/portal/data_page/) and [26Q1 release notes](https://forum.depmap.org/t/announcing-the-26q1-release/4606). The cache selects matched cell models from CRISPRGeneEffect.csv and joins model metadata. It adds gene aliases, dependency labels, and an unmatched-model row. Cite the specific DepMap release and portal. Keep the applicable [DepMap terms](https://depmap.org/portal/terms/) with any further use; this project does not assign CC BY or MIT to the 26Q1 data. The portal terms page required interactive verification during preparation.
- **Eide2023:** [publication](https://doi.org/10.1158/2643-3230.BCD-23-0014), Table S9 drug responses and [GSE216087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216087) guide counts. The app derives exploratory guide-level ratios and gene summaries; see the normalization limitation in METHODS.md.
- **Zhang2024:** the local `41419_2024_6995_MOESM2_ESM.xlsx` file supplies guide-library design context. The source publication is [CRISPR screen of venetoclax response-associated genes identifies transcription factor ZNF740 as a key functional regulator](https://doi.org/10.1038/s41419-024-06995-x). The context layer is not an independently scored screen.

## What the public snapshot contains

The repository contains project code, curated gene sets, and source/evidence metadata. Prepared app caches are distributed separately as assets on release v0.1.3. Raw source PDFs, spreadsheets, and full downloaded expression/count matrices are omitted. Numerical scores have been harmonized and in some cases reanalyzed, so they should not be described as verbatim original supplements. Source score conventions and transformation limitations are described in [METHODS.md](METHODS.md).

The [cache manifest](../data/metadata/cache_manifest.csv) records SHA-256 hashes and byte sizes of the exact distributed cache files. The [Nakao analysis manifest](../data/metadata/nakao2023_analysis_manifest.json) preserves the original source/result hashes, versions, and commands with project-root paths made relative. The Nakao intermediate gene_results.tsv is regenerated by the documented builder and is not needed to run the app.

## Third-party code

The vendored DrugZ implementation is copyright 2018 hart-lab under its [MIT license](../scripts/vendor/DRUGZ_LICENSE). MAGeCK is installed separately from its pinned upstream source archive; it is not vendored or relicensed here.

## Citation practice

Cite the studies supporting a result, not just this browser. Cite the dataset release for DepMap and the BeatAML2 manuscript/data repository when using their context layers. Candidate articles in the additional-screen catalog are literature leads; they are not all incorporated datasets.
