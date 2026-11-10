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

## v0.1.2 release (3 November 2026)

v0.1.2 is the first release after v0.1.0; a v0.1.1 cache set was prepared but never published. The caches were rebuilt from local source inputs. `depmap_crispr.rds` is unchanged from v0.1.0; `gene_scores.rds`, `external_context.rds`, and `additional_screen_candidates.rds` were rebuilt.

- Gene-score cache: 884,653 rows, 16 studies, 68 displayed experiments (v0.1.0: 439,550 rows, 10 studies, 37 experiments). Added studies: Wang2023, Sharon2019, Yokoyama2026, Stewart2026, BrimEdwards2026, and Diepstraten2026; Deng2022 now contributes its CRISPRa hit lists instead of RNA-seq differential expression.
- Fixes relative to v0.1.0: MAGeCK selection arms keep the signed gene LFC (Glytsou2023, Nechiporuk2019), Eide2023 count summaries are depth-normalized and median-centered, Portelinha2024 ranks are computed across genes, and Stewart2026's unnormalized venetoclax arms are median-centered with their published FDRs kept aside. For each change, rows outside the intended studies were checked to be byte-identical before and after.
- `scripts/smoke_test.R` passed against the source rebuild, including biological checks for every added study (for example BAX enrichment in all MAGeCK and Eide2023 venetoclax arms, A1/BCL2A1 enrichment in both Deng2022 CRISPRa lines, and no FDR-based significance from the Stewart2026 unnormalized arms).
- `Rscript --vanilla scripts/check_release.R` passed with the updated cache manifest, including Gene Explorer screen-view, metric, and expansion queries.
- BrimEdwards2026 gene results are re-derived from GEO counts by `scripts/build_brimedwards2026.R`, which needs edgeR (Bioconductor; not in `renv.lock`). The distributed cache already contains these results, so running the app does not need edgeR.

### Verification after publication (v0.1.2)

Release v0.1.2 was published on 3 November 2026 from `main` at c37e764 (merges of #1 and #2; commit IDs as rewritten on 3 November 2026 to remove attribution trailers), with the four caches above as assets. In a fresh clone of `main` with no local caches:

- `scripts/download_prepared_data.R` downloaded all four assets from v0.1.2, and each matched the byte size and SHA-256 recorded in `data/metadata/cache_manifest.csv`. A second run skipped all four as already verified.
- `Rscript --vanilla scripts/check_release.R` passed: 884,653 rows, 16 studies, 68 displayed experiments, including app construction and the Gene Explorer screen-view, metric, and expansion server queries.

Not checked: a clean dependency installation on another operating system, or a source rebuild of the BrimEdwards2026 results with edgeR.

## v0.1.3 release (10 November 2026)

v0.1.3 adds baseline DepMap context and redesigns the Paper Explorer. Only `depmap_crispr.rds` changed; `gene_scores.rds`, `external_context.rds`, and `additional_screen_candidates.rds` are byte-identical to v0.1.2.

- `depmap_crispr.rds` was rebuilt from DepMap Public 26Q1 `CRISPRGeneEffect.csv` with MOLM-14 (ACH-001574), MV4-11 (ACH-000045), and OCI-AML2 (ACH-000113) added to MOLM-13, THP-1, and OCI-LY19: 111,186 rows. Gene effects for the three original lines match v0.1.2 to within 1e-14 (CSV round-trip); one gene key changed under the current alias rules (ADSS1 is keyed ADSSL1).
- `Rscript --vanilla scripts/check_release.R` passed with the updated manifest, including Paper Explorer screen buttons for every paper and result, and the Gene Explorer baseline dependency rows (MCL1 in the DepMap tail, group off when unchecked).
- `scripts/smoke_test.R` was not rerun: it forces a full source rebuild, and only the DepMap cache changed.
