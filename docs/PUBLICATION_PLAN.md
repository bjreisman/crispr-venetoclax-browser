# GitHub release plan

Target: public repository `bjreisman/crispr-venetoclax-browser`, with MIT-licensed original code and separately attributed public datasets distributed as GitHub Release assets. This first version shares the existing exploratory project; it does not claim to resolve all scientific review findings.

## 1. Prepare the first public version

- [x] Choose repository visibility, name, and original-code license with the project owner.
- [x] Make normal startup read prepared caches without triggering source rebuilds.
- [x] Exclude raw downloads, source PDFs/spreadsheets, local credentials, virtual environments, deployment metadata, and internal review materials.
- [x] Add methods/limitations, quick-start instructions, data-source records, attribution, and software citation metadata.
- [x] Capture installed R package versions in `renv.lock`.
- [x] Generate cache byte sizes and SHA-256 hashes.
- [x] Export an explicit file list to a clean publication checkout, inspect it for secrets and oversized files, and preserve the original working folder/history.
- [x] Run cache, loader, and app-startup checks in a directory containing only public files.
- [x] Create the public repository, push the code snapshot, publish v0.1.0 with four cache assets, and verify the download script and app from a fresh GitHub clone.

The local `.git/config` is a OneDrive cloud placeholder and local Git reads time out. If it remains unavailable, use a clean initial public snapshot instead of modifying or uploading unverified historical commits. The research folder and existing Git history stay intact.

Published: [repository](https://github.com/bjreisman/crispr-venetoclax-browser) and [v0.1.0 data release](https://github.com/bjreisman/crispr-venetoclax-browser/releases/tag/v0.1.0). All four assets were downloaded from GitHub into a fresh clone and matched the recorded hashes; the app/server checks and independent localhost startup passed.

## 2. Improve reproducibility after the initial upload

- Automate acquisition of source inputs using stable download locations, versions, and source checksums; the first release documents manual acquisition.
- Check the locked dependency installation on another machine or CI runner. Recording a lockfile is separate from proving cross-platform installation.
- Record dated analysis releases and changelogs. Keep binary caches in versioned release assets, or move them to a data archive if distribution needs change.
- Add an automated check for cache/schema compatibility and startup when updating code or prepared data.

## 3. Scientific improvements to pursue separately

1. Define assay/treatment eligibility and study-level directional support using the existing screen/experiment registries. Count independent studies and models separately from tests, timepoints, and aggregates.
2. Separate adjusted significance, nominal evidence, source thresholds, missing statistics, and measured unsupported results. Make missing-FDR filtering explicit.
3. Document verified contrast orientations (MAGeCK selection-arm signs are preserved as of v0.1.2). Make heatmaps and recurrence tables use one evidence definition, with source-native effect units in detail views.
4. Reanalyze Eide counts with replicate-aware scoring and guide/control QC (depth normalization and median centering were added in v0.1.2).
5. Add species-aware gene identity and versioned orthology; keep source symbols and ambiguous mappings visible.
6. Handle repeated specimens/subjects in patient associations, and ensure empty selections and exports have explicit behavior.

These scientific tasks remain open and are disclosed in the app and methods. They are not prerequisites for sharing a clearly labeled exploratory project.
