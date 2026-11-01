# Initial public snapshot validation

Checked 1 November 2026 on macOS with R 4.6.1.

- All four distributed app caches match the byte sizes and SHA-256 values in `data/metadata/cache_manifest.csv`.
- Gene-score cache: 439,550 rows, ten studies, 37 displayed experiment labels. Labels do not represent independent screens.
- `Rscript --vanilla scripts/check_release.R` passed in the research folder and the clean public export using installed packages matching the recorded lockfile. Checks cover expected data fields, controls, DepMap models, external context, candidate catalog consistency, cache-only loading, app construction, and representative gene/overview/expression server queries.
- A regression fixture with an old cache timestamp loaded the prepared data without calling builders. Missing caches produced actionable errors without writing replacements.
- The app started from the public export with no raw source folders and served HTTP 200 with navigation and interpretation notices. The local HTTP check ran outside the sandbox because its listening sockets are blocked.
- The code export excludes binary caches as well as raw study folders, original PDFs/spreadsheets, virtual environments, deployment-account directory, internal review notes, or inherited Git history.
- Text files and character fields in the binary caches were checked for common GitHub/AWS/private-key credential patterns and local home-directory paths; none were found. This pattern check is not a guarantee that arbitrary secrets cannot exist.
- Largest release asset: `external_context.rds`, 76,132,852 bytes. All four caches are excluded from Git history and distributed as versioned release assets.

A pre-existing ggplot warning about a dropped `text` aesthetic appears during server checks; the checks pass. This preparation does not audit all source-paper contrasts, rerun the scientific analyses, or validate target-prioritization results. Restoring all dependencies from the network on another machine remains a separate reproducibility check.

## Verification after publication

The public repository and release v0.1.0 were verified on 1 November 2026. All four release assets were downloaded by `scripts/download_prepared_data.R` into a fresh clone from GitHub and passed the recorded SHA-256 and size checks. The release check script passed using the installed packages matching the lockfile. The freshly cloned app started independently on localhost port 3844 and returned HTTP 200. No `.rds` files are tracked in Git; the data are release assets.

The downloader was also checked to skip already verified local files and to preserve an existing mismatched cache while reporting an actionable error. These checks do not establish a clean dependency installation on another operating system.
