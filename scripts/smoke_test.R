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

if (!nrow(depmap_crispr)) {
  stop("DepMap CRISPR data are missing. Run scripts/build_depmap_data.R.")
}

expected_depmap_lines <- c("MOLM-13", "THP-1", "OCI-LY-19")
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
