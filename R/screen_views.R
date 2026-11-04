library(dplyr)
library(tibble)

# One representative result per independent screen for the Gene Explorer
# summary view. Repeat timepoints use the latest timepoint; MAGeCK screens list
# both selection arms under one screen_id, and each gene takes the arm that
# matches its LFC sign. Parent screens follow
# data/metadata/experiment_to_screen_mapping.csv.
screen_view_groups <- c(
  "Venetoclax · AML",
  "Venetoclax · lymphoid",
  "Venetoclax combinations",
  "Other BH3 mimetics",
  "Other assays (direction not comparable)"
)

venetoclax_screen_groups <- screen_view_groups[1:3]

screen_view_registry <- tibble::tribble(
  ~screen_id, ~paper, ~experiment, ~screen_group, ~short_label, ~hit_list_only,
  "chen_molm13_ven", "Chen2019", "Venetoclax vs DMSO day 16", "Venetoclax · AML", "Chen 2019 · MOLM-13", FALSE,
  "lin_molm13_ven", "Lin2019", "MOLM-13 drug sensitivity", "Venetoclax · AML", "Lin 2019 · MOLM-13 (metabolic library)", FALSE,
  "lin_thp1_ven", "Lin2019", "THP-1 drug sensitivity", "Venetoclax · AML", "Lin 2019 · THP-1 (metabolic library)", FALSE,
  "nechiporuk_brunello_ven", "Nechiporuk2019", "Brunello library venetoclax vs DMSO day 14 | positive selection", "Venetoclax · AML", "Nechiporuk 2019 · MOLM-13 Brunello", FALSE,
  "nechiporuk_brunello_ven", "Nechiporuk2019", "Brunello library venetoclax vs DMSO day 14 | negative selection", "Venetoclax · AML", "Nechiporuk 2019 · MOLM-13 Brunello", FALSE,
  "nechiporuk_ykosuke_ven", "Nechiporuk2019", "Y. Kosuke library venetoclax vs DMSO day 14 | positive selection", "Venetoclax · AML", "Nechiporuk 2019 · MOLM-13 Yusa", FALSE,
  "eide_ven", "Eide2023", "OCI-AML2 Yusa library day 21 venetoclax vs DMSO", "Venetoclax · AML", "Eide 2023 · OCI-AML2", FALSE,
  "stewart_ven", "Stewart2026", "OCI-AML2 Yusa library day 21 venetoclax 200 nM vs DMSO", "Venetoclax · AML", "Stewart 2026 · OCI-AML2", FALSE,
  "stewart_ven_palbo", "Stewart2026", "OCI-AML2 Yusa library day 21 venetoclax 200 nM + palbociclib 1 µM vs DMSO", "Venetoclax combinations", "Stewart 2026 · OCI-AML2 VEN + palbociclib", FALSE,
  "nakao_wt_ven", "Nakao2023", "March5-WT mMA9Cas venetoclax 3.5 µM vs DMSO", "Venetoclax · AML", "Nakao 2023 · mouse MLL-AF9", FALSE,
  "nakao_ko_ven", "Nakao2023", "March5-KO mMA9Cas venetoclax 3 µM vs DMSO", "Venetoclax · AML", "Nakao 2023 · mouse MLL-AF9, March5-KO", FALSE,
  "wang_brunello_ven", "Wang2023", "Brunello genome-wide library venetoclax 25 nM vs DMSO day 20", "Venetoclax · AML", "Wang 2023 · MOLM-13 Brunello", FALSE,
  "wang_rbp_ven", "Wang2023", "RBP domain-focused library venetoclax vs DMSO", "Venetoclax · AML", "Wang 2023 · MOLM-13 (RBP domain library)", FALSE,
  "sharon_molm13r1_ven", "Sharon2019", "MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day 29 | positive selection", "Venetoclax · AML", "Sharon 2019 · MOLM13-R1 (acquired resistant)", FALSE,
  "sharon_molm13r1_ven", "Sharon2019", "MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day 29 | negative selection", "Venetoclax · AML", "Sharon 2019 · MOLM13-R1 (acquired resistant)", FALSE,
  "guieze_ko_ven", "Guièze2019", "STARS replicate 1", "Venetoclax · lymphoid", "Guièze 2019 · OCI-Ly1 (hit list)", TRUE,
  "portelinha_ven", "Portelinha2024", "Positive selection screen table 1", "Venetoclax · lymphoid", "Portelinha 2024 · OCI-LY19 (hit list)", TRUE,
  "yokoyama_ed_ven", "Yokoyama2026", "ED40515(-) venetoclax vs DMSO | positive selection", "Venetoclax · lymphoid", "Yokoyama 2026 · PTCL ED40515(−)", FALSE,
  "yokoyama_ed_ven", "Yokoyama2026", "ED40515(-) venetoclax vs DMSO | negative selection", "Venetoclax · lymphoid", "Yokoyama 2026 · PTCL ED40515(−)", FALSE,
  "yokoyama_del_ven", "Yokoyama2026", "DEL venetoclax vs DMSO | positive selection", "Venetoclax · lymphoid", "Yokoyama 2026 · PTCL DEL", FALSE,
  "yokoyama_del_ven", "Yokoyama2026", "DEL venetoclax vs DMSO | negative selection", "Venetoclax · lymphoid", "Yokoyama 2026 · PTCL DEL", FALSE,
  "yokoyama_ed_navitoclax", "Yokoyama2026", "ED40515(-) navitoclax vs DMSO | positive selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL ED40515(−) navitoclax", FALSE,
  "yokoyama_ed_navitoclax", "Yokoyama2026", "ED40515(-) navitoclax vs DMSO | negative selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL ED40515(−) navitoclax", FALSE,
  "yokoyama_del_navitoclax", "Yokoyama2026", "DEL navitoclax vs DMSO | positive selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL DEL navitoclax", FALSE,
  "yokoyama_del_navitoclax", "Yokoyama2026", "DEL navitoclax vs DMSO | negative selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL DEL navitoclax", FALSE,
  "yokoyama_ed_abbv467", "Yokoyama2026", "ED40515(-) ABBV-467 vs DMSO | positive selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL ED40515(−) ABBV-467 (MCL1i)", FALSE,
  "yokoyama_ed_abbv467", "Yokoyama2026", "ED40515(-) ABBV-467 vs DMSO | negative selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL ED40515(−) ABBV-467 (MCL1i)", FALSE,
  "yokoyama_del_abbv467", "Yokoyama2026", "DEL ABBV-467 vs DMSO | positive selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL DEL ABBV-467 (MCL1i)", FALSE,
  "brimedwards_azd4320", "BrimEdwards2026", "OCI-AML2 Yusa library day 14 AZD4320 30 nM vs DMSO", "Other BH3 mimetics", "Brim-Edwards 2026 · OCI-AML2 AZD4320 (BCL2/BCL-XLi)", FALSE,
  "yokoyama_del_abbv467", "Yokoyama2026", "DEL ABBV-467 vs DMSO | negative selection", "Other BH3 mimetics", "Yokoyama 2026 · PTCL DEL ABBV-467 (MCL1i)", FALSE,
  "eide_ven_rux", "Eide2023", "OCI-AML2 Yusa library day 21 venetoclax+ruxolitinib vs DMSO", "Venetoclax combinations", "Eide 2023 · OCI-AML2 VEN + ruxolitinib", FALSE,
  "glytsou_ven_aza", "Glytsou2023", "V+A | positive selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + azacitidine", FALSE,
  "glytsou_ven_aza", "Glytsou2023", "V+A | negative selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + azacitidine", FALSE,
  "glytsou_mcl1i_ven", "Glytsou2023", "MCL1i+V | positive selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + MCL1i", FALSE,
  "glytsou_mcl1i_ven", "Glytsou2023", "MCL1i+V | negative selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + MCL1i", FALSE,
  "glytsou_mcl1i", "Glytsou2023", "MCL1i | positive selection", "Other BH3 mimetics", "Glytsou 2023 · MOLM-13 MCL1i", FALSE,
  "glytsou_mcl1i", "Glytsou2023", "MCL1i | negative selection", "Other BH3 mimetics", "Glytsou 2023 · MOLM-13 MCL1i", FALSE,
  "nakao_wt_s63845", "Nakao2023", "March5-WT mMA9Cas S63845 15 µM vs DMSO", "Other BH3 mimetics", "Nakao 2023 · mouse MLL-AF9 S63845", FALSE,
  "portelinha_s63845", "Portelinha2024", "Positive selection screen table 2", "Other BH3 mimetics", "Portelinha 2024 · OCI-LY19 S63845 (hit list)", TRUE,
  "deng_crispra_214", "Deng2022", "CRISPRa line #214 venetoclax IC80 (30 nM) vs DMSO | enriched sgRNAs", "Other assays (direction not comparable)", "Deng 2022 · mouse lymphoma #214 CRISPRa (hit list)", TRUE,
  "deng_crispra_216", "Deng2022", "CRISPRa line #216 venetoclax IC80 (15 nM) vs DMSO | enriched sgRNAs", "Other assays (direction not comparable)", "Deng 2022 · mouse lymphoma #216 CRISPRa (hit list)", TRUE,
  "diepstraten_crispra_214", "Diepstraten2026", "CRISPRa Partita DHL214 venetoclax IC50 (10 nM) vs DMSO | sgRNAs >=5% frequency (hit list)", "Other assays (direction not comparable)", "Diepstraten 2026 · mouse lymphoma DHL214 Partita CRISPRa (hit list)", TRUE,
  "diepstraten_crispra_214", "Diepstraten2026", "CRISPRa Partita DHL214 venetoclax IC80 (30 nM) vs DMSO | sgRNAs >=5% frequency (hit list)", "Other assays (direction not comparable)", "Diepstraten 2026 · mouse lymphoma DHL214 Partita CRISPRa (hit list)", TRUE,
  "diepstraten_crispra_216", "Diepstraten2026", "CRISPRa Partita DHL216 venetoclax IC50 (5 nM) vs DMSO | sgRNAs >=5% frequency (hit list)", "Other assays (direction not comparable)", "Diepstraten 2026 · mouse lymphoma DHL216 Partita CRISPRa (hit list)", TRUE,
  "diepstraten_crispra_216", "Diepstraten2026", "CRISPRa Partita DHL216 venetoclax IC80 (15 nM) vs DMSO | sgRNAs >=5% frequency (hit list)", "Other assays (direction not comparable)", "Diepstraten 2026 · mouse lymphoma DHL216 Partita CRISPRa (hit list)", TRUE,
  "guieze_orf_ven", "Guièze2019", "ORF overexpression screen | venetoclax 100 nM day 10", "Other assays (direction not comparable)", "Guièze 2019 · ORF overexpression", FALSE,
  "he_p1r", "He2026", "P1-R BH3 priming screen CytoC-high vs CytoC-low", "Other assays (direction not comparable)", "He 2026 · AML PDX P1-R BH3 priming", FALSE,
  "he_p2r", "He2026", "P2-R BH3 priming screen CytoC-high vs CytoC-low", "Other assays (direction not comparable)", "He 2026 · AML PDX P2-R BH3 priming", FALSE
)

# Top or bottom 10% of a screen on the signed percentile scale (-1 to +1).
screen_signal_cutoff <- 0.8

# Reported FDRs of 0 are drawn at this signed -log10 value.
screen_fdr_cap <- 50

# Axis options for the screen view. `band` shades |value| >= band; percentile
# and robust z are undefined for hit-list tables, which get their own column.
screen_metrics <- list(
  percentile = list(
    label = "Within-screen percentile",
    column = "percentile",
    band = screen_signal_cutoff,
    band_note = "Shaded: top or bottom 10% of the screen",
    axis = "← sensitizing    within-screen percentile    resistance →",
    hit_list_column = TRUE
  ),
  robust_z = list(
    label = "Robust z-score",
    column = "robust_z",
    band = 3,
    band_note = "Shaded: |robust z| ≥ 3 (median and MAD within the screen)",
    axis = "← sensitizing    robust z-score within screen    resistance →",
    hit_list_column = TRUE
  ),
  signed_fdr = list(
    label = "Signed −log₁₀ FDR",
    column = "signed_log10_fdr",
    band = -log10(0.05),
    band_note = paste0("Shaded: FDR ≤ 0.05, as published. FDR reported as 0 is drawn at ±", screen_fdr_cap, ". Stewart2026 venetoclax and combination FDRs are shown as published but not counted as significant."),
    axis = "← sensitizing    signed −log₁₀ FDR as published (directional arm)    resistance →",
    hit_list_column = FALSE
  ),
  raw = list(
    label = "Raw effect (as published)",
    column = "published_effect",
    band = NA_real_,
    band_note = "Raw effects as published, in each paper's own units; not comparable between screens. Stewart2026 venetoclax and combination values are unnormalized as published (median about −6).",
    axis = "← sensitizing    raw effect value as published (paper units)    resistance →",
    hit_list_column = FALSE
  )
)

screen_metric_choices <- stats::setNames(names(screen_metrics), vapply(screen_metrics, `[[`, character(1), "label"))

build_screen_scores <- function(scores) {
  for (col in c("published_effect_value", "published_fdr")) {
    if (!col %in% names(scores)) scores[[col]] <- NA_real_
  }
  scores |>
    inner_join(screen_view_registry, by = c("paper", "experiment")) |>
    filter(!is.na(.data$effect_value)) |>
    mutate(
      arm = dplyr::case_when(
        grepl("\\| positive selection$", .data$experiment) ~ "positive",
        grepl("\\| negative selection$", .data$experiment) ~ "negative",
        TRUE ~ "single"
      ),
      paired = .data$screen_id %in% .data$screen_id[.data$arm == "negative"]
    ) |>
    # Both MAGeCK arms carry the same LFC, so exactly one arm matches its sign.
    filter(
      !.data$paired |
        (.data$arm == "positive" & .data$effect_value >= 0) |
        (.data$arm == "negative" & .data$effect_value < 0)
    ) |>
    # Mouse paralogs merged under one human key keep their most extreme value.
    arrange(desc(abs(.data$effect_value))) |>
    distinct(.data$screen_id, .data$gene_key, .keep_all = TRUE) |>
    group_by(.data$screen_id) |>
    mutate(
      percentile = dplyr::if_else(
        .data$hit_list_only,
        NA_real_,
        (rank(.data$effect_value, ties.method = "average") - 0.5) / dplyr::n() * 2 - 1
      ),
      screen_genes = dplyr::n(),
      robust_z = dplyr::if_else(
        .data$hit_list_only | stats::mad(.data$effect_value) == 0,
        NA_real_,
        (.data$effect_value - stats::median(.data$effect_value)) / stats::mad(.data$effect_value)
      )
    ) |>
    ungroup() |>
    mutate(
      # Published values; they differ from effect_value/fdr only where an
      # importer recentred an unnormalized screen (Stewart2026).
      published_effect = dplyr::coalesce(.data$published_effect_value, .data$effect_value),
      published_fdr = dplyr::coalesce(.data$published_fdr, .data$fdr),
      signed_log10_fdr = dplyr::case_when(
        is.na(.data$published_fdr) ~ NA_real_,
        .data$published_fdr <= 0 ~ sign(.data$published_effect) * screen_fdr_cap,
        TRUE ~ sign(.data$published_effect) * pmin(-log10(.data$published_fdr), screen_fdr_cap)
      ),
      supported = (!is.na(.data$fdr) & .data$fdr <= 0.05) | grepl("^paper threshold", .data$significance),
      screen_group = factor(.data$screen_group, levels = screen_view_groups)
    ) |>
    select(
      "screen_id", "screen_group", "short_label", "hit_list_only", "paper", "experiment",
      "experiment_id", "gene_key", "gene_symbol", "effect_metric", "effect_value",
      "p_value", "fdr", "significance", "supported", "percentile", "robust_z",
      "published_effect", "published_fdr", "signed_log10_fdr", "screen_genes"
    )
}

summarise_screen_support <- function(screen_scores, gene_keys, papers = NULL) {
  screens <- screen_view_registry |>
    filter(.data$screen_group %in% venetoclax_screen_groups) |>
    distinct(.data$screen_id, .data$paper, .data$hit_list_only)
  if (!is.null(papers)) screens <- screens |> filter(.data$paper %in% papers)
  rows <- screen_scores |>
    filter(.data$screen_id %in% screens$screen_id, .data$gene_key %in% gene_keys)
  purrr::map_dfr(gene_keys, function(key) {
    gene_rows <- rows |> filter(.data$gene_key == key)
    # Hit-list screens count as tested for every gene: absence means not a hit.
    tested <- union(gene_rows$screen_id, screens$screen_id[screens$hit_list_only])
    resistance <- gene_rows |>
      filter(dplyr::coalesce(.data$percentile >= screen_signal_cutoff, .data$hit_list_only & .data$supported))
    sensitizing <- gene_rows |>
      filter(!is.na(.data$percentile), .data$percentile <= -screen_signal_cutoff)
    tibble::tibble(
      gene_key = key,
      gene_symbol = dplyr::first(c(gene_rows$gene_symbol, key)),
      n_tested = length(tested),
      n_resistance = nrow(resistance),
      n_resistance_fdr = sum(resistance$supported),
      n_sensitizing = nrow(sensitizing),
      n_sensitizing_fdr = sum(sensitizing$supported)
    )
  })
}
