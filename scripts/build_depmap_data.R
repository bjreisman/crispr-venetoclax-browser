source("R/load_data.R")

depmap_dir <- file.path("data", "external", "depmap")
selected_effect_file <- file.path(depmap_dir, "CRISPRGeneEffect_selected_26Q1.csv")
model_file <- file.path(depmap_dir, "Model_26Q1.csv")
out_file <- file.path("data", "processed", "depmap_crispr.rds")

if (!file.exists(selected_effect_file)) {
  stop("Missing ", selected_effect_file, ". Extract matched rows from CRISPRGeneEffect_26Q1.csv first.")
}

if (!file.exists(model_file)) {
  stop("Missing ", model_file, ".")
}

model_map <- readr::read_csv(model_file, show_col_types = FALSE) |>
  dplyr::transmute(
    depmap_model_id = .data$ModelID,
    cell_line = .data$CellLineName,
    stripped_cell_line = .data$StrippedCellLineName,
    ccle_name = .data$CCLEName,
    oncotree_lineage = .data$OncotreeLineage,
    oncotree_primary_disease = .data$OncotreePrimaryDisease
  )

study_cell_lines <- tibble::tribble(
  ~study_cell_line, ~depmap_model_id, ~matching_note,
  "MOLM-13", "ACH-000362", "direct match",
  "THP-1", "ACH-000146", "direct match",
  "OCI-LY19", "ACH-000124", "direct match",
  "OCI-Ly1", NA_character_, "not present in DepMap Public 26Q1 Model.csv"
)

effect_matrix <- readr::read_csv(selected_effect_file, show_col_types = FALSE)
names(effect_matrix)[1] <- "depmap_model_id"

depmap_crispr <- effect_matrix |>
  tidyr::pivot_longer(
    cols = -depmap_model_id,
    names_to = "gene",
    values_to = "gene_effect"
  ) |>
  dplyr::mutate(
    gene_symbol = stringr::str_remove(.data$gene, " \\([0-9]+\\)$"),
    entrez_id = stringr::str_match(.data$gene, "\\(([0-9]+)\\)$")[, 2],
    gene_key = gene_key(.data$gene_symbol),
    depmap_release = "DepMap Public 26Q1",
    source_dataset = "CRISPRGeneEffect.csv",
    dependency_strength = dplyr::case_when(
      is.na(.data$gene_effect) ~ "not measured",
      .data$gene_effect <= -1 ~ "strong dependency",
      .data$gene_effect <= -0.5 ~ "moderate dependency",
      TRUE ~ "weak/no dependency"
    )
  ) |>
  dplyr::left_join(model_map, by = "depmap_model_id") |>
  dplyr::left_join(study_cell_lines, by = "depmap_model_id") |>
  dplyr::select(
    depmap_release, source_dataset, study_cell_line, matching_note,
    depmap_model_id, cell_line, stripped_cell_line, ccle_name,
    oncotree_lineage, oncotree_primary_disease,
    gene_symbol, gene_key, entrez_id, gene_effect, dependency_strength
  )

unmatched <- study_cell_lines |>
  dplyr::filter(is.na(.data$depmap_model_id)) |>
  dplyr::mutate(
    depmap_release = "DepMap Public 26Q1",
    source_dataset = "CRISPRGeneEffect.csv",
    depmap_model_id = NA_character_,
    cell_line = .data$study_cell_line,
    stripped_cell_line = NA_character_,
    ccle_name = NA_character_,
    oncotree_lineage = NA_character_,
    oncotree_primary_disease = NA_character_,
    gene_symbol = NA_character_,
    gene_key = NA_character_,
    entrez_id = NA_character_,
    gene_effect = NA_real_,
    dependency_strength = "not measured"
  ) |>
  dplyr::select(dplyr::all_of(names(depmap_crispr)))

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(dplyr::bind_rows(depmap_crispr, unmatched), out_file)

message("Built ", out_file)
message("Rows: ", nrow(depmap_crispr))
message("Matched DepMap models: ", paste(sort(unique(depmap_crispr$cell_line)), collapse = ", "))
