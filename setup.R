required <- c("shiny", "bslib", "tidyverse", "readxl", "DT", "plotly", "ggrepel")
optional <- c("shinyWidgets", "renv")

missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
missing_optional <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_required)) {
  install.packages(missing_required, repos = "https://cloud.r-project.org")
}

message("Required package check complete.")
if (length(missing_optional)) {
  message("Optional packages not installed: ", paste(missing_optional, collapse = ", "))
  message("Install them with install.packages(c(", paste(sprintf('\"%s\"', missing_optional), collapse = ", "), ")) if you want renv or shinyWidgets later.")
}
