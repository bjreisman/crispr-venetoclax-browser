source("R/load_data.R")

nakao_source <- "Nakao2023/41375_2023_1879_MOESM2_ESM.xlsx"
nakao_result <- "data/processed/nakao2023/gene_results.tsv"
nakao_script <- "scripts/build_nakao2023_analysis.py"
nakao_python <- Sys.getenv("NAKAO2023_PYTHON", unset = ".venv-nakao2023/bin/python")

nakao_needs_build <- file.exists(nakao_source) &&
  (!file.exists(nakao_result) || file.info(nakao_result)$mtime < file.info(nakao_source)$mtime)

if (nakao_needs_build) {
  if (!file.exists(nakao_python)) {
    stop(
      "Nakao2023 analysis environment is missing. Create .venv-nakao2023, install scripts/requirements-nakao2023.txt, then rerun this script."
    )
  }
  status <- system2(nakao_python, c(nakao_script, "--root", "."))
  if (!identical(status, 0L)) stop("Nakao2023 MAGeCK/DrugZ analysis failed; browser data were not rebuilt.")
}

scores <- load_browser_data(force = TRUE)
external_context <- write_external_context(force = TRUE)

message("Built data/processed/gene_scores.rds and data/processed/gene_scores.csv")
message("Rows: ", nrow(scores))
message("Papers: ", paste(sort(unique(scores$paper)), collapse = ", "))
message("Experiments: ", length(unique(scores$experiment_id)))
message("Built data/processed/external_context.rds")
message("Drug-response rows: ", nrow(external_context$drug_response))
message("BeatAML2 expression genes: ", nrow(external_context$beataml_expression$genes))
message("BeatAML2 expression samples: ", ncol(external_context$beataml_expression$expression))
message("BeatAML2 gene sets: ", nrow(external_context$beataml_gene_sets$coverage))
message("Eide CRISPR context rows: ", nrow(external_context$eide_crispr))
message("Zhang library rows: ", nrow(external_context$zhang_library))
