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
  "nakao_wt_ven", "Nakao2023", "March5-WT mMA9Cas venetoclax 3.5 µM vs DMSO", "Venetoclax · AML", "Nakao 2023 · mouse MLL-AF9", FALSE,
  "nakao_ko_ven", "Nakao2023", "March5-KO mMA9Cas venetoclax 3 µM vs DMSO", "Venetoclax · AML", "Nakao 2023 · mouse MLL-AF9, March5-KO", FALSE,
  "guieze_ko_ven", "Guièze2019", "STARS replicate 1", "Venetoclax · lymphoid", "Guièze 2019 · OCI-Ly1 (hit list)", TRUE,
  "portelinha_ven", "Portelinha2024", "Positive selection screen table 1", "Venetoclax · lymphoid", "Portelinha 2024 · OCI-LY19 (hit list)", TRUE,
  "eide_ven_rux", "Eide2023", "OCI-AML2 Yusa library day 21 venetoclax+ruxolitinib vs DMSO", "Venetoclax combinations", "Eide 2023 · OCI-AML2 VEN + ruxolitinib", FALSE,
  "glytsou_ven_aza", "Glytsou2023", "V+A | positive selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + azacitidine", FALSE,
  "glytsou_ven_aza", "Glytsou2023", "V+A | negative selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + azacitidine", FALSE,
  "glytsou_mcl1i_ven", "Glytsou2023", "MCL1i+V | positive selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + MCL1i", FALSE,
  "glytsou_mcl1i_ven", "Glytsou2023", "MCL1i+V | negative selection", "Venetoclax combinations", "Glytsou 2023 · MOLM-13 VEN + MCL1i", FALSE,
  "glytsou_mcl1i", "Glytsou2023", "MCL1i | positive selection", "Other BH3 mimetics", "Glytsou 2023 · MOLM-13 MCL1i", FALSE,
  "glytsou_mcl1i", "Glytsou2023", "MCL1i | negative selection", "Other BH3 mimetics", "Glytsou 2023 · MOLM-13 MCL1i", FALSE,
  "nakao_wt_s63845", "Nakao2023", "March5-WT mMA9Cas S63845 15 µM vs DMSO", "Other BH3 mimetics", "Nakao 2023 · mouse MLL-AF9 S63845", FALSE,
  "portelinha_s63845", "Portelinha2024", "Positive selection screen table 2", "Other BH3 mimetics", "Portelinha 2024 · OCI-LY19 S63845 (hit list)", TRUE,
  "guieze_orf_ven", "Guièze2019", "ORF overexpression screen | venetoclax 100 nM day 10", "Other assays (direction not comparable)", "Guièze 2019 · ORF overexpression", FALSE,
  "he_p1r", "He2026", "P1-R BH3 priming screen CytoC-high vs CytoC-low", "Other assays (direction not comparable)", "He 2026 · AML PDX P1-R BH3 priming", FALSE,
  "he_p2r", "He2026", "P2-R BH3 priming screen CytoC-high vs CytoC-low", "Other assays (direction not comparable)", "He 2026 · AML PDX P2-R BH3 priming", FALSE,
  "deng_expression", "Deng2022", "Eµ-Myc/dCas9a-SAM vs B-ALL expression", "Other assays (direction not comparable)", "Deng 2022 · mouse lymphoma expression", FALSE
)

# Top or bottom 10% of a screen on the signed percentile scale (-1 to +1).
screen_signal_cutoff <- 0.8

build_screen_scores <- function(scores) {
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
    distinct(.data$screen_id, .data$gene_key, .keep_all = TRUE) |>
    group_by(.data$screen_id) |>
    mutate(
      percentile = dplyr::if_else(
        .data$hit_list_only,
        NA_real_,
        (rank(.data$effect_value, ties.method = "average") - 0.5) / dplyr::n() * 2 - 1
      ),
      screen_genes = dplyr::n()
    ) |>
    ungroup() |>
    mutate(
      supported = (!is.na(.data$fdr) & .data$fdr <= 0.05) | grepl("^paper threshold", .data$significance),
      screen_group = factor(.data$screen_group, levels = screen_view_groups)
    ) |>
    select(
      "screen_id", "screen_group", "short_label", "hit_list_only", "paper", "experiment",
      "experiment_id", "gene_key", "gene_symbol", "effect_metric", "effect_value",
      "p_value", "fdr", "significance", "supported", "percentile", "screen_genes"
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
