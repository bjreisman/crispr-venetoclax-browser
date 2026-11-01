library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

trapezoid_auc <- function(concentration, response) {
  keep <- is.finite(concentration) & is.finite(response) & concentration > 0
  concentration <- concentration[keep]
  response <- response[keep]
  if (length(unique(concentration)) < 2) return(NA_real_)
  x <- log10(concentration)
  ord <- order(x)
  x <- x[ord]
  y <- response[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / diff(range(x))
}

tidy_beataml2_drug_response <- function(root = ".") {
  source_file <- file.path(root, "data/external/beataml2/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt")
  if (!file.exists(source_file)) return(tibble())

  readr::read_tsv(source_file, show_col_types = FALSE) |>
    filter(.data$passed_qc, .data$inhibitor == "Venetoclax") |>
    group_by(
      dataset = "BeatAML2",
      sample_id = as.character(.data$dbgap_subject_id),
      rnaseq_sample = .data$dbgap_rnaseq_sample,
      dnaseq_sample = .data$dbgap_dnaseq_sample,
      drug = .data$inhibitor,
      treatment_context = .data$inhibitor_panel,
      time_of_read = .data$time_of_read
    ) |>
    summarise(
      response_metric = "computed viability AUC from raw dose response",
      response_value = trapezoid_auc(.data$well_concentration, .data$normalized_viability),
      min_viability = min(.data$normalized_viability, na.rm = TRUE),
      n_doses = n_distinct(.data$well_concentration),
      source_file = source_file,
      .groups = "drop"
    ) |>
    mutate(
      response_direction = "lower AUC = more sensitive",
      disease_context = "Primary AML specimens"
    )
}

tidy_beataml2_expression <- function(root = ".") {
  source_file <- file.path(root, "data/external/beataml2/beataml_waves1to4_norm_exp_dbgap.txt")
  if (!file.exists(source_file)) {
    return(list(
      expression = matrix(numeric(), nrow = 0, ncol = 0),
      genes = tibble(gene_key = character(), gene_symbol = character(), description = character(), biotype = character())
    ))
  }

  raw <- readr::read_tsv(source_file, show_col_types = FALSE)
  sample_cols <- setdiff(names(raw), c("stable_id", "display_label", "description", "biotype"))

  expression <- raw |>
    mutate(
      gene_symbol = normalize_gene_symbol(.data$display_label),
      gene_key = gene_key(.data$gene_symbol)
    ) |>
    filter(!is.na(.data$gene_key), .data$biotype == "protein_coding") |>
    group_by(.data$gene_key) |>
    summarise(
      gene_symbol = dplyr::first(stats::na.omit(.data$gene_symbol)),
      description = dplyr::first(stats::na.omit(.data$description)),
      biotype = dplyr::first(stats::na.omit(.data$biotype)),
      across(all_of(sample_cols), ~ mean(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  genes <- expression |> select("gene_key", "gene_symbol", "description", "biotype")
  mat <- expression |>
    select(all_of(sample_cols)) |>
    as.matrix()
  storage.mode(mat) <- "numeric"
  rownames(mat) <- expression$gene_key

  list(expression = mat, genes = genes)
}

load_curated_gene_sets <- function(root = ".") {
  source_file <- file.path(root, "data/gene_sets_curated.csv")
  if (!file.exists(source_file)) {
    return(tibble(gene_set = character(), category = character(), gene_symbol = character(), gene_key = character()))
  }

  readr::read_csv(source_file, show_col_types = FALSE) |>
    transmute(
      gene_set = as.character(.data$gene_set),
      category = as.character(.data$category),
      gene_symbol = normalize_gene_symbol(.data$gene_symbol),
      gene_key = gene_key(.data$gene_symbol)
    ) |>
    filter(!is.na(.data$gene_set), !is.na(.data$gene_key)) |>
    distinct(.data$gene_set, .data$category, .data$gene_symbol, .data$gene_key)
}

score_beataml_gene_sets <- function(beataml_expression, gene_sets) {
  expression_matrix <- beataml_expression$expression
  empty <- list(
    definitions = gene_sets,
    scores = matrix(numeric(), nrow = 0, ncol = 0),
    coverage = tibble(
      gene_set = character(),
      category = character(),
      n_genes_defined = integer(),
      n_genes_scored = integer(),
      genes_scored = character(),
      genes_missing = character()
    )
  )
  if (!length(expression_matrix) || !nrow(gene_sets)) return(empty)

  present_genes <- intersect(unique(gene_sets$gene_key), rownames(expression_matrix))
  if (!length(present_genes)) return(empty)

  scaled_expression <- t(scale(t(expression_matrix[present_genes, , drop = FALSE])))
  scaled_expression[!is.finite(scaled_expression)] <- NA_real_

  coverage <- gene_sets |>
    group_by(.data$gene_set, .data$category) |>
    summarise(
      n_genes_defined = n_distinct(.data$gene_key),
      n_genes_scored = n_distinct(.data$gene_key[.data$gene_key %in% rownames(scaled_expression)]),
      genes_scored = paste(sort(unique(.data$gene_symbol[.data$gene_key %in% rownames(scaled_expression)])), collapse = ", "),
      genes_missing = paste(sort(unique(.data$gene_symbol[!.data$gene_key %in% rownames(scaled_expression)])), collapse = ", "),
      .groups = "drop"
    ) |>
    arrange(.data$category, .data$gene_set)

  scored_sets <- coverage |> filter(.data$n_genes_scored > 0)
  if (!nrow(scored_sets)) return(empty)

  score_rows <- lapply(scored_sets$gene_set, function(set_name) {
    genes <- gene_sets |>
      filter(.data$gene_set == set_name, .data$gene_key %in% rownames(scaled_expression)) |>
      pull(.data$gene_key) |>
      unique()
    colMeans(scaled_expression[genes, , drop = FALSE], na.rm = TRUE)
  })
  scores <- do.call(rbind, score_rows)
  rownames(scores) <- scored_sets$gene_set
  storage.mode(scores) <- "numeric"

  list(definitions = gene_sets, scores = scores, coverage = coverage)
}

tidy_eide2023_drug_response <- function(root = ".") {
  source_file <- file.path(root, "data/external/eide2023/bcd-23-0014_table_s9_suppst9.xlsx")
  if (!file.exists(source_file)) return(tibble())

  suppressWarnings(readxl::read_excel(source_file, sheet = "Supp Table S9")) |>
    filter(stringr::str_detect(.data$drug, "^VEN( |$|\\+)")) |>
    transmute(
      dataset = "Eide2023",
      sample_id = as.character(.data$patient_id),
      lab_id = .data$lab_id,
      beataml_dbgap_subject_id = as.character(.data$beataml_dbgap_subject_id),
      rnaseq_sample = .data$beataml_dbgap_rnaseq_sample,
      drug = .data$drug,
      treatment_context = .data$drug_type,
      response_metric = "published ex vivo AUC (% maximum)",
      response_value = .data$auc_pct_of_max,
      probit_auc = .data$probit_auc,
      response_direction = "lower AUC = more sensitive",
      disease_context = .data$dxAtSpecimenAcquisition,
      disease_stage = .data$diseaseStageAtSpecimenCollection,
      specimen_type = .data$specimenType,
      eln_2022 = .data$`2022ELN`,
      variant_summary = .data$variantSummary,
      percent_blasts_bm = .data$percentBlastsBM,
      percent_monocytes_pb = .data$percentMonocytesPB,
      source_file = source_file
    )
}

tidy_eide2023_crispr_context <- function(root = ".") {
  scores <- tidy_eide2023_crispr(root)
  if (!nrow(scores)) return(tibble())
  scores |>
    select(
      dataset = paper,
      experiment,
      gene_symbol,
      gene_key,
      effect_metric,
      effect_value,
      rank,
      direction,
      source_file,
      any_of("n_sgrna")
    )
}

tidy_zhang2024_library <- function(root = ".") {
  source_file <- file.path(root, "data/external/zhang2024/41419_2024_6995_MOESM2_ESM.xlsx")
  if (!file.exists(source_file)) return(tibble())

  suppressWarnings(readxl::read_excel(source_file, sheet = "Table S1", skip = 2, col_names = c("sgrna_id", "sgrna_sequence", "gene_symbol"))) |>
    filter(!is.na(.data$gene_symbol), .data$gene_symbol != "gene_symbol") |>
    mutate(
      dataset = "Zhang2024",
      gene_symbol = normalize_gene_symbol(.data$gene_symbol),
      gene_key = gene_key(.data$gene_symbol),
      source_file = source_file
    ) |>
    group_by(.data$dataset, .data$gene_symbol, .data$gene_key, .data$source_file) |>
    summarise(
      n_sgrna = n(),
      example_sgrnas = paste(head(.data$sgrna_id, 4), collapse = ", "),
      .groups = "drop"
    ) |>
    arrange(.data$gene_symbol)
}

build_external_context <- function(root = ".") {
  beataml_expression <- tidy_beataml2_expression(root)
  gene_sets <- load_curated_gene_sets(root)
  list(
    drug_response = bind_rows(
      tidy_beataml2_drug_response(root),
      tidy_eide2023_drug_response(root)
    ),
    beataml_expression = beataml_expression,
    beataml_gene_sets = score_beataml_gene_sets(beataml_expression, gene_sets),
    eide_crispr = tidy_eide2023_crispr_context(root),
    zhang_library = tidy_zhang2024_library(root)
  )
}
