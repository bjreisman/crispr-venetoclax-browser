# Download the prepared data for the initial release. Does not rebuild analyses.
# Run from the repository root after installing/restoring the R dependencies.
download_prepared_data <- function(overwrite = FALSE) {
  manifest <- readr::read_csv("data/metadata/cache_manifest.csv", show_col_types = FALSE)
  base_url <- "https://github.com/bjreisman/crispr-venetoclax-browser/releases/download/v0.1.3"
  expected <- file.path("data/processed", c(
    "gene_scores.rds", "depmap_crispr.rds", "external_context.rds", "additional_screen_candidates.rds"
  ))
  stopifnot(setequal(manifest$path, expected), !anyDuplicated(manifest$path))
  valid <- function(path, row) {
    file.exists(path) && file.info(path)$size == row$bytes &&
      identical(digest::digest(file = path, algo = "sha256"), row$sha256)
  }
  fetch <- function(row) {
    target <- row$path
    if (valid(target, row)) {
      message("Already verified: ", target)
      return(invisible(NULL))
    }
    if (file.exists(target) && !overwrite) {
      stop("Existing cache differs from this release: ", target,
           ". Back it up, then use download_prepared_data(overwrite = TRUE) to replace it.")
    }
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("cache-download-", tmpdir = dirname(target))
    on.exit(unlink(temporary), add = TRUE)
    status <- utils::download.file(paste0(base_url, "/", basename(target)), temporary, mode = "wb", method = "libcurl")
    if (status != 0L || !valid(temporary, row)) {
      stop("Download failed or checksum did not match for ", target, ". Existing files were preserved.")
    }
    if (!file.copy(temporary, target, overwrite = overwrite)) stop("Could not save verified cache: ", target)
    message("Downloaded and verified: ", target)
  }
  for (i in seq_len(nrow(manifest))) fetch(manifest[i, ])
  invisible(manifest)
}

download_prepared_data()
