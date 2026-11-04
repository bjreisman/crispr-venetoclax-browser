library(rsconnect)

app_name <- Sys.getenv("SHINYAPPS_APP_NAME", unset = "crispr-venetoclax-browser")

if (!file.exists("data/processed/gene_scores.rds")) {
  stop("Missing data/processed/gene_scores.rds. Run scripts/build_processed_data.R first.")
}

if (!file.exists("data/processed/depmap_crispr.rds")) {
  stop("Missing data/processed/depmap_crispr.rds. Run scripts/build_depmap_data.R first.")
}

if (!file.exists("data/processed/external_context.rds")) {
  stop("Missing data/processed/external_context.rds. Run scripts/build_processed_data.R first.")
}

app_files <- c(
  "app.R",
  "README.md",
  "R/load_data.R",
  "R/plots.R",
  "R/screen_views.R",
  "data/metadata/experiment_to_screen_mapping.csv",
  "R/tidy_external_context.R",
  "R/tidy_sources.R",
  "data/processed/gene_scores.rds",
  "data/processed/depmap_crispr.rds",
  "data/processed/external_context.rds",
  "data/processed/additional_screen_candidates.rds"
)

missing <- app_files[!file.exists(app_files)]
if (length(missing)) {
  stop("Missing deployment files: ", paste(missing, collapse = ", "))
}

rsconnect::deployApp(
  appDir = ".",
  appName = app_name,
  appTitle = "CRISPR Venetoclax Resistance Browser",
  appFiles = app_files,
  forceUpdate = TRUE
)
