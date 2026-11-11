# CRISPR Venetoclax Browser

An R Shiny project for exploring published venetoclax resistance and sensitivity datasets across hematologic malignancy models. Browse genes, experiments, cross-study coverage, baseline dependency, and expression/drug-response context.

**Status: exploratory work in progress.** The browser combines different assays and score scales. Its recurrence counts represent dataset presence, not validated replicated sensitization. See [methods and known limitations](docs/METHODS.md) before using it to prioritize targets.

## Run with the prepared data

Clone or download this repository, open R in the project directory, then run:

```r
install.packages("renv") # once, if needed
renv::restore(prompt = FALSE)
renv::load()
source("scripts/download_prepared_data.R") # once: downloads and verifies the release caches
shiny::runApp(".")
```

The dependency snapshot was generated with R 4.6.1. Run `renv::load()` again in each new R session to use the project library. System libraries required by R packages must be installed for your platform. The four prepared `.rds` caches are distributed with [release v0.1.4](https://github.com/bjreisman/crispr-venetoclax-browser/releases/tag/v0.1.4). The download script checks byte sizes and SHA-256 hashes against the versioned manifest, and skips already verified files. No raw supplements, Python environment, or analysis rebuild are needed. Startup reads the downloaded caches without rebuilding them based on timestamps.

For an unpinned installation, `source("setup.R")` installs missing required packages from CRAN. Use `renv::restore()` for the recorded versions.

From a shell with R on PATH:

```sh
Rscript --vanilla -e 'renv::load(); source("scripts/check_release.R")'
Rscript --vanilla -e 'renv::load(); shiny::runApp(".", host="127.0.0.1", port=3838)'
```

GitHub hosts the source and data; running the interactive app requires R/Shiny or a separate Shiny hosting service.

## What is included

| Path | Purpose |
| --- | --- |
| `app.R` | Shiny UI and server |
| `R/` | Data loading, source importers, plotting |
| `data/processed/*.rds` | Four app caches downloaded from the release, approximately 87 MiB combined; excluded from Git |
| `data/metadata/` | Screen registry, experiment mapping, candidate catalog, source and cache manifests |
| `data/gene_sets_curated.csv` | Curated expression-context gene sets |
| `scripts/` | Explicit data builds, release checks, export, and deployment |
| `renv.lock` | R dependency snapshot |
| `docs/` | Methods, data attribution, and release plan |

The repository excludes prepared binary caches, raw paper folders, PDFs, source spreadsheets, large raw downloads, virtual environments, deployment account metadata, and internal review notes. Source-detail views need the corresponding local supplementary files; without them the app displays an availability note.

The approximately 73 MiB external-context cache is a GitHub Release asset, along with the other three caches. Binary data are kept out of Git history. The script preserves existing files that differ from the release unless replacement is explicitly requested.

## Rebuild deliberately from source

Rebuilding is optional and requires the original inputs at the paths in [the source manifest](data/metadata/source_manifest.csv). The initial public version is runnable from caches; obtaining all raw inputs is currently a manual step.

For the Nakao2023 reanalysis, create an isolated Python environment:

```sh
uv venv .venv-nakao2023 --python 3.11
uv pip install --python .venv-nakao2023/bin/python -r scripts/requirements-nakao2023.txt
.venv-nakao2023/bin/python scripts/build_nakao2023_analysis.py --root .
```

Then, in R from the repository root:

```r
source("scripts/build_processed_data.R")
source("scripts/build_depmap_data.R")
source("scripts/build_additional_screen_catalog.R")
```

The DepMap builder expects `Model_26Q1.csv` and either `CRISPRGeneEffect_selected_26Q1.csv` (the matched model rows) or the full `CRISPRGeneEffect_26Q1.csv` under `data/external/depmap/`; it extracts the matched rows from the full matrix when a study model is missing. It does not download DepMap files automatically.

`scripts/smoke_test.R` performs the existing full analysis checks and **forces data rebuilding**. It requires raw inputs and the Nakao environment. `scripts/check_release.R` checks prepared caches and app startup without rebuilding or changing scientific results.

## Deploy

The app is published on [Posit Connect Cloud](https://bjreisman-crispr-venetoclax-browser.share.connect.posit.cloud). Install `rsconnect` (1.6 or later), log in once with `rsconnect::connectCloudUser()`, and run:

```r
source("scripts/deploy.R")
```

Set `DEPLOY_APP_ID` to update an existing content item and `DEPLOY_APP_NAME` to choose a different name. For shinyapps.io, configure credentials with `rsconnect::setAccountInfo()` and run `source("scripts/deploy_shinyapps.R")`. Keep tokens and secrets out of repository files. The deployment bundle includes app code and all four prepared caches: gene scores, DepMap, external context, and additional-screen candidates. It excludes raw source files.

## Data, license, and citation

Original project code is provided under the [MIT License](LICENSE). Third-party data retain their upstream terms; the MIT license does not relicense datasets. See [DATA_SOURCES.md](docs/DATA_SOURCES.md) for attribution, source links, and transformations. The vendored DrugZ implementation retains its [original MIT notice](scripts/vendor/DRUGZ_LICENSE).

For use in research, cite the underlying studies and dataset releases as well as this software. Repository citation metadata are provided in `CITATION.cff`.

The preparation checklist and future scientific improvements are tracked in [PUBLICATION_PLAN.md](docs/PUBLICATION_PLAN.md).
