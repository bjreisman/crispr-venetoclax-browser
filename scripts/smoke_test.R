source("R/load_data.R")

scores <- load_browser_data(force = TRUE)
depmap_crispr <- load_depmap_crispr()
external_context <- load_external_context(force = TRUE)

required <- standard_columns
missing <- setdiff(required, names(scores))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

paper_counts <- table(scores$paper)
if (any(paper_counts == 0) || length(paper_counts) < 10) {
  stop("Expected nonempty rows for all core papers plus Eide2023, He2026, and Nakao2023. Counts: ", paste(names(paper_counts), paper_counts, sep = "=", collapse = ", "))
}

known_hits <- c("BAX", "TP53", "PMAIP1", "CLPB", "MFN2", "SETD1B", "MCL1", "BCL2", "ADSS2")
present_hits <- known_hits[gene_key(known_hits) %in% scores$gene_key]
message("Known hits present: ", paste(present_hits, collapse = ", "))

if (!all(c("BAX", "TP53", "PMAIP1") %in% present_hits)) {
  stop("Core positive-control hits BAX, TP53, and PMAIP1 were not all found.")
}

nakao <- scores |> dplyr::filter(.data$paper == "Nakao2023")
expected_nakao_experiments <- c(
  "March5-WT mMA9Cas venetoclax 3.5 µM vs DMSO",
  "March5-WT mMA9Cas S63845 15 µM vs DMSO",
  "March5-KO mMA9Cas venetoclax 3 µM vs DMSO"
)
if (!all(expected_nakao_experiments %in% nakao$experiment) ||
    dplyr::n_distinct(nakao$experiment) != length(expected_nakao_experiments)) {
  stop("Nakao2023 should contain exactly three drug-vs-DMSO experiments.")
}

if (any(grepl("^(NonTargetingControlGuideForMouse_|BRDN|mmu-mir-)", nakao$gene_symbol, ignore.case = TRUE))) {
  stop("Nakao2023 browser rows should exclude control and miRNA pseudo-genes.")
}

if (any(!is.finite(nakao$p_value) | nakao$p_value < 0 | nakao$p_value > 1) ||
    any(!is.finite(nakao$fdr) | nakao$fdr < 0 | nakao$fdr > 1)) {
  stop("Nakao2023 DrugZ p-values and browser-facing FDRs must lie in [0, 1].")
}

nakao_effect <- function(experiment_value, genes) {
  nakao |>
    dplyr::filter(.data$experiment == experiment_value, .data$gene_key %in% gene_key(genes)) |>
    dplyr::select(.data$gene_key, .data$effect_value, .data$mageck_beta_delta)
}

wt_ven <- nakao_effect(expected_nakao_experiments[1], c("MARCH5", "UBE2J2", "UBE2K"))
if (nrow(wt_ven) != 3 || any(wt_ven$effect_value >= 0) || any(wt_ven$mageck_beta_delta >= 0)) {
  stop("Nakao2023 March5-WT venetoclax screen should deplete MARCH5, UBE2J2, and UBE2K by DrugZ and MAGeCK.")
}

wt_s63845 <- nakao_effect(expected_nakao_experiments[2], c("BAX", "BAK1", "APAF1", "CASP9"))
if (nrow(wt_s63845) != 4 || any(wt_s63845$effect_value <= 0) || any(wt_s63845$mageck_beta_delta <= 0)) {
  stop("Nakao2023 S63845 screen should enrich apoptosis-pathway knockouts by DrugZ and MAGeCK.")
}

ko_ven <- nakao_effect(expected_nakao_experiments[3], c("PMAIP1", "BAX", "APAF1", "CASP9", "TRP53"))
if (nrow(ko_ven) != 5 || any(ko_ven$effect_value <= 0) || any(ko_ven$mageck_beta_delta <= 0)) {
  stop("Nakao2023 March5-KO venetoclax screen should enrich PMAIP1, BAX, APAF1, CASP9, and TP53 knockouts.")
}

ko_e2 <- nakao_effect(expected_nakao_experiments[3], c("UBE2J2", "UBE2K"))
if (nrow(ko_e2) != 2 || any(ko_e2$effect_value < -1)) {
  stop("Nakao2023 March5-KO screen should not retain strong UBE2J2/UBE2K synthetic-lethal effects.")
}

nakao_detail <- load_detail_for_gene("MARCH5", "Nakao2023")
if (!nrow(nakao_detail) || !all(c("sgrna_id", "source_sheet") %in% names(nakao_detail))) {
  stop("Nakao2023 source detail should return raw sgRNA rows from the supplement.")
}

adss_family <- scores |>
  dplyr::filter(.data$gene_symbol %in% c("ADSS", "Adss", "ADSS2", "ADSSL1", "Adssl1"))

if (!all(adss_family$gene_key[adss_family$gene_symbol %in% c("ADSS", "Adss", "ADSS2")] == "ADSS2")) {
  stop("ADSS/ADSS2 aliases are not aligned to ADSS2.")
}

if (!all(adss_family$gene_key[adss_family$gene_symbol %in% c("ADSSL1", "Adssl1")] == "ADSSL1")) {
  stop("ADSSL1 should remain separate from ADSS2.")
}

if (!all(scores$gene_key[scores$gene_symbol %in% c("TOMM70", "TOMM70A", "Tomm70a")] == "TOMM70")) {
  stop("TOMM70/TOMM70A aliases are not aligned.")
}

if (!all(gene_key(c("Bcl2a1a", "Bcl2a1b", "Bcl2a1c", "Bcl2a1d", "BCL2A1")) == "BCL2A1")) {
  stop("Mouse A1 paralogs should map to human BCL2A1.")
}

deng_crispra <- scores |>
  dplyr::filter(.data$paper == "Deng2022", grepl("^CRISPRa line #21[46] ", .data$experiment))

if (dplyr::n_distinct(deng_crispra$experiment) != 2L || anyDuplicated(deng_crispra[c("experiment", "gene_key")])) {
  stop("Deng2022 CRISPRa should contain two enriched-sgRNA hit lists with one row per gene key.")
}

deng_a1 <- deng_crispra |> dplyr::filter(.data$gene_key == "BCL2A1")
if (nrow(deng_a1) != 2L || any(deng_a1$fdr > 0.05) || any(deng_a1$effect_value <= 0)) {
  stop("Deng2022 CRISPRa should enrich A1 (BCL2A1) at FDR <= 0.05 in both lines.")
}

if (!all(c("BCL2L1", "MCL1") %in% deng_crispra$gene_key[grepl("#216", deng_crispra$experiment)])) {
  stop("Deng2022 CRISPRa line #216 should enrich BCL2L1 and MCL1.")
}

orf_hits <- scores |>
  dplyr::filter(
    .data$paper == "Guièze2019",
    grepl("ORF overexpression", .data$experiment),
    .data$gene_key %in% gene_key(c("BCL2L1", "BCL2L2", "BCL2", "MCL1"))
  )

if (nrow(orf_hits) != 4) {
  stop("Guièze2019 ORF screen did not include expected BCL2-family top hits.")
}

orf_threshold_hits <- scores |>
  dplyr::filter(
    .data$paper == "Guièze2019",
    grepl("ORF overexpression", .data$experiment),
    .data$significance == "paper threshold: ORF LFC >= 2"
  )

if (nrow(orf_threshold_hits) != 70) {
  stop("Guièze2019 ORF paper-threshold gene count should match the reported 70 genes.")
}

lin_threshold_hits <- scores |>
  dplyr::filter(
    .data$paper == "Lin2019",
    .data$experiment == "Average drug sensitivity",
    .data$significance == "paper threshold: Avg_sens < -0.914"
  )

if (nrow(lin_threshold_hits) != 241) {
  stop("Lin2019 paper-threshold average drug-sensitivity hit count should match the reported 241 genes.")
}

direction_sign_mismatches <- scores |>
  dplyr::filter(
    (.data$direction == "resistance/enriched" & .data$effect_value < 0) |
      (.data$direction == "sensitization/depleted" & .data$effect_value > 0)
  )

if (nrow(direction_sign_mismatches)) {
  stop("Direction labels disagree with effect signs in ", nrow(direction_sign_mismatches), " rows.")
}

# BAX loss confers venetoclax resistance; it must be enriched in both MAGeCK
# selection arms and in the Eide2023 count-derived screens.
bax_not_enriched <- scores |>
  dplyr::filter(
    .data$gene_key == "BAX",
    grepl("^MAGeCK (positive|negative) selection LFC$", .data$effect_metric) | .data$paper == "Eide2023",
    grepl("V|venetoclax", .data$experiment),
    !(.data$effect_value > 0)
  )

if (nrow(bax_not_enriched)) {
  stop("BAX is not enriched in: ", paste(unique(bax_not_enriched$experiment_id), collapse = "; "))
}

eide_offsets <- scores |>
  dplyr::filter(.data$paper == "Eide2023") |>
  dplyr::group_by(.data$experiment_id) |>
  dplyr::summarise(median_effect = stats::median(.data$effect_value, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(abs(.data$median_effect) > 0.1)

if (nrow(eide_offsets)) {
  stop("Eide2023 fold-changes are not centered: ", paste(eide_offsets$experiment_id, collapse = "; "))
}

constant_ranks <- scores |>
  dplyr::filter(!is.na(.data$rank)) |>
  dplyr::group_by(.data$experiment_id) |>
  dplyr::summarise(n = dplyr::n(), n_ranks = dplyr::n_distinct(.data$rank), .groups = "drop") |>
  dplyr::filter(.data$n > 1, .data$n_ranks == 1)

if (nrow(constant_ranks)) {
  stop("Experiments with a constant rank: ", paste(constant_ranks$experiment_id, collapse = "; "))
}

# Wang2023 smoke checks; append to scripts/smoke_test.R (operates on `scores`).
wang <- scores |> dplyr::filter(.data$paper == "Wang2023")
expected_wang_experiments <- c(
  "Brunello genome-wide library venetoclax 25 nM vs DMSO day 20",
  "RBP domain-focused library venetoclax vs DMSO"
)
if (!setequal(unique(wang$experiment), expected_wang_experiments) ||
    sum(wang$experiment == expected_wang_experiments[1]) < 19000 ||
    sum(wang$experiment == expected_wang_experiments[2]) < 450) {
  stop("Wang2023 should contain the genome-wide (~19,115 genes) and RBP domain (~491 genes) venetoclax arms only.")
}

if (any(!is.na(wang$p_value)) || any(!is.na(wang$fdr))) {
  stop("Wang2023 Table S1 publishes no p-values or FDRs; they must stay NA.")
}

wang_gw <- wang |> dplyr::filter(.data$experiment == expected_wang_experiments[1])
wang_gw_top <- wang_gw |> dplyr::slice_max(.data$effect_value, n = 3)
if (!setequal(wang_gw_top$gene_key, c("PMAIP1", "BAX", "TP53")) ||
    wang_gw$effect_value[wang_gw$gene_key == "BCL2L1"] >= 0) {
  stop("Wang2023 genome-wide venetoclax arm should rank PMAIP1, BAX and TP53 as the top enriched knockouts and deplete BCL2L1.")
}

wang_rbp <- wang |>
  dplyr::filter(.data$experiment == expected_wang_experiments[2], .data$gene_key %in% c("RBM10", "SRSF11", "SRSF8", "HNRNPD", "HNRNPAB", "HNRNPF"))
if (nrow(wang_rbp) != 6 || any(wang_rbp$effect_value >= 0)) {
  stop("Wang2023 RBP domain screen should deplete the validated venetoclax sensitizers RBM10, SRSF11, SRSF8, HNRNPD, HNRNPAB and HNRNPF.")
}

sharon <- scores |> dplyr::filter(.data$paper == "Sharon2019")
sharon_days <- c(8, 16, 22, 29)
expected_sharon_experiments <- as.vector(outer(
  paste0("MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day ", sharon_days),
  c(" | negative selection", " | positive selection"),
  paste0
))
if (!setequal(unique(sharon$experiment), expected_sharon_experiments)) {
  stop("Sharon2019 should contain both MAGeCK arms for days 8, 16, 22 and 29 only.")
}

if (any(sharon$gene_symbol %in% c("chr10Promi", "chr10Rand", "EGFP", "LacZ", "luciferase"))) {
  stop("Sharon2019 browser rows should exclude TKO control targets.")
}

sharon_effect <- function(day, arm, genes) {
  sharon |>
    dplyr::filter(
      .data$experiment == paste0("MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day ", day, " | ", arm, " selection"),
      .data$gene_key %in% gene_key(genes)
    )
}

bax <- sharon_effect(29, "positive", "BAX")
if (nrow(bax) != 1 || bax$effect_value <= 0 || bax$rank != 1 || bax$fdr > 0.05) {
  stop("Sharon2019 day 29 positive selection should rank BAX first with positive LFC and FDR <= 0.05.")
}

# Fig. 1B: four of the top-10 day-8 negatively selected genes are mitochondrial translation genes.
mito_d8 <- sharon_effect(8, "negative", c("DAP3", "MRPL54", "MRPL17", "RBFA"))
if (nrow(mito_d8) != 4 || any(mito_d8$effect_value >= 0) || any(mito_d8$rank > 10)) {
  stop("Sharon2019 day 8 negative selection should rank DAP3, MRPL54, MRPL17 and RBFA in the top 10 with negative LFC.")
}

mito_d29 <- sharon |>
  dplyr::filter(grepl("day 29 \\| negative", .data$experiment), grepl("^MRP[LS][0-9]+$", .data$gene_key))
if (nrow(mito_d29) < 50 || stats::median(mito_d29$effect_value) > -0.5) {
  stop("Sharon2019 day 29 mitoribosome knockouts should be depleted overall (median LFC <= -0.5).")
}

yokoyama <- scores |> dplyr::filter(.data$paper == "Yokoyama2026")

if (dplyr::n_distinct(yokoyama$experiment) != 12L || anyDuplicated(yokoyama[c("experiment", "gene_key")])) {
  stop("Yokoyama2026 should contain six model-drug groups x two MAGeCK arms with one row per gene key.")
}

if (any(!is.na(yokoyama$fdr)) || any(yokoyama$p_value < 0 | yokoyama$p_value > 1, na.rm = TRUE)) {
  stop("Yokoyama2026 has MAGeCK p-values in [0, 1] and no published FDRs.")
}

if (any(grepl("^Control_|^[0-9]{1,2}-(Mar|Sep|Dec)$", yokoyama$gene_symbol)) ||
    !all(c("MARCH5", "SEPT9", "SEP15") %in% yokoyama$gene_key)) {
  stop("Yokoyama2026 should drop non-targeting controls and restore Excel-dated symbols from RefSeq IDs.")
}

yokoyama_effect <- function(experiment_value, genes) {
  yokoyama |>
    dplyr::filter(.data$experiment == paste(experiment_value, "| positive selection"), .data$gene_key %in% gene_key(genes))
}

navitoclax_mylip <- dplyr::bind_rows(
  yokoyama_effect("ED40515(-) navitoclax vs DMSO", "MYLIP"),
  yokoyama_effect("DEL navitoclax vs DMSO", "MYLIP")
)
if (nrow(navitoclax_mylip) != 2L || any(navitoclax_mylip$effect_value >= 0)) {
  stop("Yokoyama2026 navitoclax screens should deplete MYLIP knockouts in both PTCL lines.")
}

abbv467 <- dplyr::bind_rows(
  yokoyama_effect("ED40515(-) ABBV-467 vs DMSO", c("BAX", "BCL2L1")),
  yokoyama_effect("DEL ABBV-467 vs DMSO", c("BAX", "BCL2L1"))
)
if (nrow(abbv467) != 4L ||
    any(abbv467$effect_value[abbv467$gene_key == "BAX"] <= 0) ||
    any(abbv467$effect_value[abbv467$gene_key == "BCL2L1"] >= 0)) {
  stop("Yokoyama2026 ABBV-467 screens should enrich BAX and deplete BCL2L1 knockouts.")
}

stewart <- scores |> dplyr::filter(.data$paper == "Stewart2026")
expected_stewart_experiments <- c(
  "OCI-AML2 Yusa library day 21 venetoclax 200 nM vs DMSO",
  "OCI-AML2 Yusa library day 21 venetoclax 200 nM + palbociclib 1 µM vs DMSO",
  "OCI-AML2 Yusa library day 21 palbociclib 1 µM vs DMSO (not a BH3 mimetic)"
)
if (!setequal(unique(stewart$experiment), expected_stewart_experiments) ||
    any(dplyr::count(stewart, .data$experiment)$n < 17000)) {
  stop("Stewart2026 should contain three genome-wide (>17,000 gene) drug-vs-DMSO experiments.")
}

stewart_effect <- function(experiment_value, genes) {
  stewart |>
    dplyr::filter(.data$experiment == experiment_value, .data$gene_key %in% gene_key(genes)) |>
    dplyr::pull(.data$effect_value, name = .data$gene_key)
}

# The venetoclax and combination arms are median-centered because the published
# counts were unnormalized; their published values are kept but their p/FDR must
# not drive significance calls.
stewart_unnormalized <- stewart |> dplyr::filter(.data$experiment %in% expected_stewart_experiments[1:2])
if (any(!is.na(stewart_unnormalized$fdr)) || any(!is.na(stewart_unnormalized$p_value)) ||
    any(stewart_unnormalized$significance != "not annotated/significant") ||
    sum(stewart_unnormalized$published_fdr <= 0.05, na.rm = TRUE) < 10000) {
  stop("Stewart2026 venetoclax and combination arms should keep published FDRs aside without counting them as significant.")
}
stewart_medians <- stewart_unnormalized |>
  dplyr::group_by(.data$experiment) |>
  dplyr::summarise(median = stats::median(.data$effect_value, na.rm = TRUE), published = stats::median(.data$published_effect_value, na.rm = TRUE))
if (any(abs(stewart_medians$median) > 0.01) || any(stewart_medians$published > -4)) {
  stop("Stewart2026 venetoclax and combination arms should be median-centered, with unnormalized published values retained.")
}
stewart_ven_median <- 0
ven_resist <- stewart_effect(expected_stewart_experiments[1], c("BAX", "PMAIP1"))
ven_sens <- stewart_effect(expected_stewart_experiments[1], c("BCL2L1", "MCL1"))
if (length(ven_resist) != 2 || any(ven_resist <= 2) ||
    length(ven_sens) != 2 || any(ven_sens >= stewart_ven_median)) {
  stop("Stewart2026 venetoclax screen should enrich BAX/PMAIP1 and deplete BCL2L1/MCL1 below the screen median.")
}

combo_bax <- stewart_effect(expected_stewart_experiments[2], "BAX")
palbo_rb1 <- stewart_effect(expected_stewart_experiments[3], "RB1")
if (length(combo_bax) != 1 || combo_bax >= ven_resist[["BAX"]] - 4 ||
    length(palbo_rb1) != 1 || palbo_rb1 <= 2) {
  stop("Stewart2026: BAX enrichment should be lost under venetoclax+palbociclib and RB1 should be enriched under palbociclib.")
}

# Diepstraten2026 (Partita CRISPRa, mouse DHL lymphoma, venetoclax hit lists).
diepstraten <- scores |> filter(paper == "Diepstraten2026")
stopifnot(nrow(diepstraten) > 0)

# Four hit lists: two lines x two doses; all rows are enriched (positive) only.
stopifnot(dplyr::n_distinct(diepstraten$experiment) == 4)
stopifnot(all(diepstraten$effect_value > 0))
stopifnot(all(diepstraten$direction == "resistance/enriched"))

# No invented significance: the venetoclax screen was not statistically tested.
stopifnot(all(is.na(diepstraten$p_value)), all(is.na(diepstraten$fdr)))

# Prosurvival BCL-2 family activation is enriched, as reported in the paper.
stopifnot("BCL2L1" %in% diepstraten$gene_key, "MCL1" %in% diepstraten$gene_key)
stopifnot(any(diepstraten$gene_key == "BCL2"))

# Mouse symbols are mapped to human keys and one row per gene per experiment.
stopifnot(!any(duplicated(diepstraten[c("experiment", "gene_key")])))

brimedwards <- scores |> dplyr::filter(.data$paper == "BrimEdwards2026")
if (!identical(unique(brimedwards$experiment), "OCI-AML2 Yusa library day 14 AZD4320 30 nM vs DMSO") ||
    nrow(brimedwards) < 17000 || anyDuplicated(brimedwards$gene_key)) {
  stop("BrimEdwards2026 should contain one genome-wide (>17,000 gene) AZD4320-vs-DMSO experiment with one row per gene.")
}

brimedwards_effect <- function(genes, column = "effect_value") {
  brimedwards |>
    dplyr::filter(.data$gene_key %in% gene_key(genes)) |>
    dplyr::pull(dplyr::all_of(column), name = .data$gene_key)
}
be_resist <- brimedwards_effect(c("BAX", "PMAIP1", "TP53", "RBM15", "ZMYND8"))
be_resist_fdr <- brimedwards_effect(c("BAX", "PMAIP1", "TP53", "RBM15", "ZMYND8"), "fdr")
be_bcl2l1 <- brimedwards_effect("BCL2L1")
if (length(be_resist) != 5 || any(be_resist <= 2) || any(be_resist_fdr > 0.05) ||
    length(be_bcl2l1) != 1 || be_bcl2l1 >= -1) {
  stop("BrimEdwards2026 AZD4320 screen should enrich BAX/PMAIP1/TP53/RBM15/ZMYND8 (FDR <= 0.05) and deplete BCL2L1.")
}

# Figure 2A labels nine Tier1 genes; the re-derived tiers should recover them.
be_tier1 <- sort(brimedwards$gene_key[brimedwards$tier %in% "Tier1"])
if (!identical(be_tier1, sort(gene_key(c("BAX", "PMAIP1", "TP53", "RBM15", "ZMYND8", "LETM2", "MRPL16", "ATL2", "PCGF1"))))) {
  stop("BrimEdwards2026 Tier1 genes should match the nine labeled in Figure 2A.")
}

if (any(brimedwards$p_value < 0 | brimedwards$p_value > 1 | brimedwards$fdr < brimedwards$p_value, na.rm = TRUE) ||
    mean(brimedwards$fdr <= 0.05, na.rm = TRUE) > 0.05) {
  stop("BrimEdwards2026 p-values/FDRs should be valid and flag only a small minority of genes.")
}

if (!nrow(depmap_crispr)) {
  stop("DepMap CRISPR data are missing. Run scripts/build_depmap_data.R.")
}

expected_depmap_lines <- c("MOLM-13", "MOLM-14", "MV4-11", "OCI-AML2", "THP-1", "OCI-LY-19")
missing_depmap_lines <- setdiff(expected_depmap_lines, unique(depmap_crispr$cell_line))
if (length(missing_depmap_lines)) {
  stop("Missing expected DepMap cell lines: ", paste(missing_depmap_lines, collapse = ", "))
}

if (!all(gene_key(c("BAX", "MCL1", "BCL2")) %in% depmap_crispr$gene_key)) {
  stop("Expected DepMap dependency genes are missing.")
}

if (!nrow(external_context$drug_response)) {
  stop("External drug-response context is missing.")
}

if (!all(c("BeatAML2", "Eide2023") %in% unique(external_context$drug_response$dataset))) {
  stop("Expected BeatAML2 and Eide2023 drug-response context.")
}

if (!nrow(external_context$beataml_expression$genes) ||
    !all(gene_key(c("BCL2", "MCL1", "BCL2L11")) %in% external_context$beataml_expression$genes$gene_key)) {
  stop("Expected BeatAML2 expression genes are missing.")
}

if (is.null(external_context$beataml_gene_sets$scores) ||
    !nrow(external_context$beataml_gene_sets$scores) ||
    !all(c("BCL2 family anti-apoptotic", "Oxidative phosphorylation", "Monocytic differentiation") %in% rownames(external_context$beataml_gene_sets$scores))) {
  stop("Expected BeatAML2 gene-set scores are missing.")
}

if (!nrow(external_context$eide_crispr) || !"Eide2023" %in% scores$paper) {
  stop("Eide2023 CRISPR context was not incorporated.")
}

he_hits <- scores |>
  dplyr::filter(
    .data$paper == "He2026",
    grepl("combined", .data$experiment),
    .data$gene_key %in% gene_key(c("ADSS2", "GART", "PFAS", "PAICS", "PPAT"))
  )

if (nrow(he_hits) != 5 || any(he_hits$effect_value >= 0)) {
  stop("He2026 BH3-priming screen did not capture expected CytoC-low enriched purine synthesis hits.")
}

he_tp53 <- scores |>
  dplyr::filter(.data$paper == "He2026", grepl("combined", .data$experiment), .data$gene_key == "TP53")

if (nrow(he_tp53) != 1 || he_tp53$effect_value <= 0) {
  stop("He2026 TP53 should be represented as CytoC-high enriched.")
}

if (any(is.na(he_hits$p_value)) || any(!is.na(he_hits$fdr))) {
  stop("He2026 should include MAGeCK p-values and leave FDR blank.")
}

he_rpm_hits <- scores |>
  dplyr::filter(
    .data$paper == "He2026",
    grepl("^P[12]-R BH3 priming screen", .data$experiment),
    .data$gene_key %in% gene_key(c("ADSS2", "GART", "PFAS", "PAICS", "PPAT"))
  )

if (nrow(he_rpm_hits) != 10 || any(he_rpm_hits$effect_value >= 0)) {
  stop("He2026 P1-R/P2-R RPM summaries did not capture expected CytoC-low enriched purine synthesis hits.")
}

if (any(!is.na(he_rpm_hits$p_value)) || any(!is.na(he_rpm_hits$fdr))) {
  stop("He2026 RPM-derived sample summaries should not include p-values or FDR.")
}

if (!"ZNF740" %in% external_context$zhang_library$gene_key) {
  stop("Zhang2024 CRISPR library context is missing ZNF740.")
}

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required.")
}

message("Smoke test passed.")
