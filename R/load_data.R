library(readr)
library(dplyr)

source("R/tidy_sources.R")
source("R/tidy_external_context.R")

processed_dir <- file.path("data", "processed")
gene_scores_rds <- file.path(processed_dir, "gene_scores.rds")
gene_scores_csv <- file.path(processed_dir, "gene_scores.csv")
depmap_crispr_rds <- file.path(processed_dir, "depmap_crispr.rds")
external_context_rds <- file.path(processed_dir, "external_context.rds")

gene_aliases <- tibble::tribble(
  ~gene_key, ~alias,
  "BCL2L11", "BIM",
  "PMAIP1", "NOXA",
  "BBC3", "PUMA",
  "MCL1", "MCL-1",
  "BCL2L1", "BCL-XL",
  "BCL2L1", "BCLX",
  "BCL2L1", "BCLXL",
  "BAK1", "BAK",
  "BID", "tBID",
  "BCL2", "BCL-2",
  "ADSS2", "ADSS",
  "ADSSL1", "ADSS1",
  "TOMM70", "TOMM70A"
)

format_gene_symbols <- function(x) {
  symbols <- sort(unique(stats::na.omit(as.character(x))))
  symbols <- symbols[order(!grepl("^[A-Z0-9.-]+$", symbols), toupper(symbols), symbols)]
  paste(symbols, collapse = " / ")
}

source_files <- c(
  "Lin2019/1-s2.0-S1550413119300117-mmc3.csv",
  "Chen2019/21598290cd190117-sup-216243_2_supp_5484296_pqktn1.xlsx",
  "Nechiporuk2019/21598290cd190125-sup-216310_2_table_5485035_pqlflj.xlsx",
  "Guieze2019/1-s2.0-S1535610819303733-mmc3.xlsx",
  "Guieze2019/1-s2.0-S1535610819303733-mmc4.xlsx",
  "Glytsou2023/cd-22-0601_supplementary_table_s1_suppst1.xlsx",
  "data/external/deng2022/41467_2022_32485_MOESM5_ESM.xlsx",
  "Portelinha2024/jem_20231143_tables1.xlsx",
  "Portelinha2024/jem_20231143_tables2.xlsx",
  "data/external/eide2023/GSE216087_crispr_counts_for_geo.txt.gz",
  "He2026/43018_2026_1184_MOESM2_ESM.xlsx",
  "He2026/43018_2026_1184_MOESM4_ESM.xlsx",
  "Nakao2023/41375_2023_1879_MOESM2_ESM.xlsx",
  "data/processed/nakao2023/gene_results.tsv",
  "data/processed/nakao2023/manifest.json"
)

cache_is_current <- function(cache_file = gene_scores_rds) {
  if (!file.exists(cache_file)) return(FALSE)
  # Invalidate when source data OR the importer code changes
  code_files <- c("R/tidy_sources.R")
  watched <- c(source_files, code_files)
  existing <- watched[file.exists(watched)]
  if (!length(existing)) return(TRUE)
  all(file.info(cache_file)$mtime >= file.info(existing)$mtime)
}

write_processed_data <- function(root = ".", force = FALSE) {
  dir.create(file.path(root, processed_dir), recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)

  if (!force && cache_is_current(gene_scores_rds)) {
    return(readRDS(gene_scores_rds))
  }

  scores <- build_gene_scores(".")
  saveRDS(scores, gene_scores_rds)
  readr::write_csv(scores, gene_scores_csv)
  scores
}

# Columns the app displays or computes with. Importers keep extra audit fields
# (DrugZ/MAGeCK side statistics, guide counts, provenance); those stay in the
# released cache but are dropped in the running app to keep memory down.
browser_score_columns <- c(
  standard_columns, "experiment_id", "significance", "signed_log10_fdr", "rationale",
  "n_orf_constructs", "n_sgrna", "published_effect_value", "published_p_value", "published_fdr"
)

slim_browser_scores <- function(scores) {
  scores[, intersect(browser_score_columns, names(scores)), drop = FALSE]
}

load_browser_data <- function(root = ".", force = FALSE) {
  if (force) return(write_processed_data(root = root, force = TRUE))
  cache <- file.path(root, gene_scores_rds)
  if (!file.exists(cache)) {
    stop("Missing ", cache, ". Restore the prepared data from the repository, or run scripts/build_processed_data.R with the source datasets available.")
  }
  readRDS(cache)
}

load_depmap_crispr <- function(root = ".") {
  depmap_file <- file.path(root, depmap_crispr_rds)
  if (!file.exists(depmap_file)) {
    return(tibble::tibble(
      depmap_release = character(),
      source_dataset = character(),
      study_cell_line = character(),
      matching_note = character(),
      depmap_model_id = character(),
      cell_line = character(),
      stripped_cell_line = character(),
      ccle_name = character(),
      oncotree_lineage = character(),
      oncotree_primary_disease = character(),
      gene_symbol = character(),
      gene_key = character(),
      entrez_id = character(),
      gene_effect = numeric(),
      dependency_strength = character()
    ))
  }
  readRDS(depmap_file)
}

write_external_context <- function(root = ".", force = FALSE) {
  dir.create(file.path(root, processed_dir), recursive = TRUE, showWarnings = FALSE)
  context_file <- file.path(root, external_context_rds)
  watched <- c(
    "R/tidy_external_context.R",
    "data/gene_sets_curated.csv",
    "data/external/beataml2/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt",
    "data/external/beataml2/beataml_waves1to4_norm_exp_dbgap.txt",
    "data/external/eide2023/bcd-23-0014_table_s9_suppst9.xlsx",
    "data/external/eide2023/GSE216087_crispr_counts_for_geo.txt.gz",
    "data/external/zhang2024/41419_2024_6995_MOESM2_ESM.xlsx"
  )
  watched <- file.path(root, watched)
  existing <- watched[file.exists(watched)]
  current <- file.exists(context_file) &&
    (force == FALSE) &&
    (!length(existing) || all(file.info(context_file)$mtime >= file.info(existing)$mtime))
  if (current) return(readRDS(context_file))
  context <- build_external_context(root)
  saveRDS(context, context_file)
  context
}

load_external_context <- function(root = ".", force = FALSE) {
  context_file <- file.path(root, external_context_rds)
  if (force) return(write_external_context(root = root, force = TRUE))
  if (!file.exists(context_file)) {
    stop("Missing ", context_file, ". Restore the prepared data from the repository, or run scripts/build_processed_data.R with the source datasets available.")
  }
  readRDS(context_file)
}

gene_choices <- function(scores) {
  choices <- scores |>
    filter(!is.na(.data$gene_key), !is.na(.data$gene_symbol)) |>
    group_by(.data$gene_key) |>
    summarise(
      gene_symbol = format_gene_symbols(.data$gene_symbol),
      papers = paste(sort(unique(.data$paper)), collapse = ", "),
      .groups = "drop"
    ) |>
    left_join(
      gene_aliases |>
        group_by(.data$gene_key) |>
        summarise(alias_label = paste(sort(unique(.data$alias)), collapse = ", "), .groups = "drop"),
      by = "gene_key"
    ) |>
    arrange(gene_symbol) |>
    mutate(
      label = if_else(
        .data$gene_key == "NONTARGETINGCONTROL",
        paste0("NTC / Non-targeting control [", .data$gene_key, "]"),
        if_else(
          is.na(.data$alias_label),
          paste0(.data$gene_symbol, " [", .data$gene_key, "]"),
          paste0(.data$gene_symbol, " / ", .data$alias_label, " [", .data$gene_key, "]")
        )
      )
    )
  stats::setNames(choices$gene_key, choices$label)
}
