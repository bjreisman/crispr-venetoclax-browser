#!/usr/bin/env python3
"""Re-analyse the Nakao2023 CRISPR screens from author-supplied sgRNA counts.

The browser displays DrugZ normZ for each drug-versus-DMSO comparison.  MAGeCK
MLE is run in parallel and retained as an audit/concordance layer; it is not
used as the browser effect score.  Outputs are intentionally compact and are
safe to commit, whereas the source workbook remains unchanged.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


WORKBOOK = Path("Nakao2023/41375_2023_1879_MOESM2_ESM.xlsx")
OUTPUT_DIR = Path("data/processed/nakao2023")
VENDORED_DRUGZ = Path("scripts/vendor/drugz.py")
EXPECTED_SHEETS = {
    "Supplemental Table 1": {"rows": 125682, "genes": 22360, "columns": 6},
    "Supplemental Table 2": {"rows": 79633, "genes": 20674, "columns": 8},
}
CONTROL_PATTERNS = (
    re.compile(r"^NonTargetingControlGuideForMouse_"),
    re.compile(r"^BRDN"),
    re.compile(r"^mmu-mir-", re.IGNORECASE),
)

ANALYSES = (
    {
        "id": "march5_wt_venetoclax",
        "sheet": "Supplemental Table 1",
        "experiment": "March5-WT mMA9Cas venetoclax 3.5 µM vs DMSO",
        "control": ["WT_output"],
        "treatment": ["WT_VEN"],
        "mle_conditions": ["WT_output", "WT_VEN", "WT_S63845"],
    },
    {
        "id": "march5_wt_s63845",
        "sheet": "Supplemental Table 1",
        "experiment": "March5-WT mMA9Cas S63845 15 µM vs DMSO",
        "control": ["WT_output"],
        "treatment": ["WT_S63845"],
        "mle_conditions": ["WT_output", "WT_VEN", "WT_S63845"],
    },
    {
        "id": "march5_ko_venetoclax",
        "sheet": "Supplemental Table 2",
        "experiment": "March5-KO mMA9Cas venetoclax 3 µM vs DMSO",
        "control": ["March5KO_output_rep1", "March5KO_output_rep2"],
        "treatment": ["March5KO_VEN_rep1", "March5KO_VEN_rep2"],
        "mle_conditions": ["March5KO_output", "March5KO_VEN"],
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def command_output(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(command, cwd=cwd, check=True, text=True, capture_output=True)
    return (completed.stdout or completed.stderr).strip()


def read_source_sheet(workbook: Path, sheet: str) -> pd.DataFrame:
    data = pd.read_excel(workbook, sheet_name=sheet, header=1, engine="openpyxl")
    expected = EXPECTED_SHEETS[sheet]
    if data.shape[0] != expected["rows"] or data.shape[1] != expected["columns"]:
        raise ValueError(f"{sheet}: expected {expected['rows']} x {expected['columns']}, found {data.shape}")
    if data["Gene"].nunique() != expected["genes"]:
        raise ValueError(f"{sheet}: expected {expected['genes']} gene labels, found {data['Gene'].nunique()}")
    if data["sgRNA"].duplicated().any() or data["sgRNA"].isna().any() or data["Gene"].isna().any():
        raise ValueError(f"{sheet}: sgRNA/Gene identifier contract failed")
    count_columns = [column for column in data.columns if column not in {"sgRNA", "Gene"}]
    if data[count_columns].isna().any().any():
        raise ValueError(f"{sheet}: count table has missing values")
    data[count_columns] = data[count_columns].astype(int)
    return data


def visible_gene(gene: str) -> bool:
    return not any(pattern.search(gene) for pattern in CONTROL_PATTERNS)


def write_count_table(data: pd.DataFrame, target: Path) -> None:
    data.to_csv(target, sep="\t", index=False, quoting=csv.QUOTE_NONE)


def write_design_matrix(spec: dict, source_columns: list[str], target: Path) -> None:
    rows: list[dict[str, int | str]] = []
    if spec["sheet"] == "Supplemental Table 1":
        condition_map = {
            "WT_input": None,
            "WT_output": "WT_output",
            "WT_VEN": "WT_VEN",
            "WT_S63845": "WT_S63845",
        }
    else:
        condition_map = {
            "March5KO_input_rep1": None,
            "March5KO_input_rep2": None,
            "March5KO_output_rep1": "March5KO_output",
            "March5KO_output_rep2": "March5KO_output",
            "March5KO_VEN_rep1": "March5KO_VEN",
            "March5KO_VEN_rep2": "March5KO_VEN",
        }
    conditions = spec["mle_conditions"]
    for sample in source_columns:
        row: dict[str, int | str] = {"Samples": sample, "baseline": 1}
        row.update({condition: 0 for condition in conditions})
        if condition_map[sample] is not None:
            row[condition_map[sample]] = 1
        rows.append(row)
    pd.DataFrame(rows).to_csv(target, sep="\t", index=False)


def run_drugz(root: Path, count_table: Path, spec: dict, raw_dir: Path) -> tuple[Path, Path, list[str]]:
    output = raw_dir / f"{spec['id']}.drugz.tsv"
    foldchange = raw_dir / f"{spec['id']}.drugz_foldchange.tsv"
    command = [
        sys.executable,
        str(VENDORED_DRUGZ),
        "-i", str(count_table),
        "-o", str(output),
        "-f", str(foldchange),
        "-c", ",".join(spec["control"]),
        "-x", ",".join(spec["treatment"]),
        "-I", "0",
    ]
    run(command, root)
    return output, foldchange, command


def run_mageck(root: Path, count_table: Path, design: Path, spec: dict, raw_dir: Path) -> tuple[Path, list[str]]:
    prefix = raw_dir / f"{spec['id']}.mageck"
    command = [
        sys.executable, "-c", "from mageck.argsParser import crisprseq_parseargs; crisprseq_parseargs()",
        "mle", "-k", str(count_table), "-d", str(design),
        "-n", str(prefix), "--norm-method", "median", "--threads", "4",
    ]
    run(command, root)
    return Path(f"{prefix}.gene_summary.txt"), command


def get_column(data: pd.DataFrame, name: str) -> pd.Series:
    return data[name] if name in data.columns else pd.Series(pd.NA, index=data.index, dtype="Float64")


def mageck_fields(mageck: pd.DataFrame, dmso_label: str, treatment_label: str) -> pd.DataFrame:
    gene_column = next((column for column in mageck.columns if column.lower() in {"id", "gene", "geneid"}), mageck.columns[0])
    result = pd.DataFrame({"Gene": mageck[gene_column].astype(str)})
    for label, prefix in ((dmso_label, "mageck_dmso"), (treatment_label, "mageck_treatment")):
        for source, target in (("beta", "beta"), ("p-value", "p_value"), ("fdr", "fdr"), ("wald-p-value", "wald_p_value"), ("wald-fdr", "wald_fdr")):
            result[f"{prefix}_{target}"] = pd.to_numeric(get_column(mageck, f"{label}|{source}"), errors="coerce")
    result["mageck_beta_delta"] = result["mageck_treatment_beta"] - result["mageck_dmso_beta"]
    return result


def browser_summary(spec: dict, source: pd.DataFrame, drugz_path: Path, mageck_path: Path) -> pd.DataFrame:
    drugz = pd.read_csv(drugz_path, sep="\t")
    if "GENE" not in drugz or "normZ" not in drugz:
        raise ValueError(f"DrugZ output schema changed: {drugz.columns.tolist()}")
    result = drugz.rename(columns={"GENE": "Gene", "numObs": "drugz_num_obs", "normZ": "effect_value"}).copy()
    result["Gene"] = result["Gene"].astype(str)
    result["effect_value"] = pd.to_numeric(result["effect_value"], errors="coerce")
    result["p_value"] = 1.0
    result["fdr"] = 1.0
    result["rank"] = pd.NA
    negative = result["effect_value"] < 0
    positive = result["effect_value"] > 0
    for direction, prefix in ((negative, "synth"), (positive, "supp")):
        result.loc[direction, "p_value"] = pd.to_numeric(result.loc[direction, f"pval_{prefix}"], errors="coerce")
        result.loc[direction, "fdr"] = pd.to_numeric(result.loc[direction, f"fdr_{prefix}"], errors="coerce").clip(upper=1)
        result.loc[direction, "rank"] = pd.to_numeric(result.loc[direction, f"rank_{prefix}"], errors="coerce")
    result["direction"] = "neutral"
    result.loc[negative, "direction"] = "sensitization/depleted"
    result.loc[positive, "direction"] = "resistance/enriched"

    guides = source.groupby("Gene", as_index=False).agg(n_sgrna=("sgRNA", "nunique"))
    mageck = pd.read_csv(mageck_path, sep="\t")
    dmso_label = "WT_output" if spec["sheet"] == "Supplemental Table 1" else "March5KO_output"
    treatment_label = spec["treatment"][0].replace("_rep1", "").replace("_rep2", "")
    result = result.merge(guides, on="Gene", how="left").merge(
        mageck_fields(mageck, dmso_label, treatment_label), on="Gene", how="left"
    )
    result = result.loc[result["Gene"].map(visible_gene)].copy()
    result["experiment"] = spec["experiment"]
    result["analysis_id"] = spec["id"]
    result["source_sheet"] = spec["sheet"]
    result["drugz_p_value_synth"] = pd.to_numeric(result["pval_synth"], errors="coerce")
    result["drugz_fdr_synth"] = pd.to_numeric(result["fdr_synth"], errors="coerce")
    result["drugz_rank_synth"] = pd.to_numeric(result["rank_synth"], errors="coerce")
    result["drugz_p_value_supp"] = pd.to_numeric(result["pval_supp"], errors="coerce")
    result["drugz_fdr_supp"] = pd.to_numeric(result["fdr_supp"], errors="coerce")
    result["drugz_rank_supp"] = pd.to_numeric(result["rank_supp"], errors="coerce")
    keep = [
        "experiment", "analysis_id", "source_sheet", "Gene", "effect_value", "p_value", "fdr", "rank", "direction",
        "n_sgrna", "drugz_num_obs", "drugz_p_value_synth", "drugz_fdr_synth", "drugz_rank_synth",
        "drugz_p_value_supp", "drugz_fdr_supp", "drugz_rank_supp", "mageck_dmso_beta", "mageck_treatment_beta",
        "mageck_beta_delta", "mageck_dmso_p_value", "mageck_dmso_fdr", "mageck_dmso_wald_p_value",
        "mageck_dmso_wald_fdr", "mageck_treatment_p_value", "mageck_treatment_fdr",
        "mageck_treatment_wald_p_value", "mageck_treatment_wald_fdr",
    ]
    return result[keep].sort_values(["experiment", "rank", "Gene"], na_position="last")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Project root")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    workbook = root / WORKBOOK
    output_dir = root / OUTPUT_DIR
    raw_dir = output_dir / "raw"
    vendor = root / VENDORED_DRUGZ
    if not workbook.exists():
        raise FileNotFoundError(workbook)
    if not vendor.exists():
        raise FileNotFoundError(vendor)
    try:
        mageck_version = command_output(
            [sys.executable, "-c", "from mageck.version import __version__; print(__version__)"], root
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError("MAGeCK is unavailable; create .venv-nakao2023 and install scripts/requirements-nakao2023.txt") from error
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    sheets = {sheet: read_source_sheet(workbook, sheet) for sheet in EXPECTED_SHEETS}
    counts: dict[str, Path] = {}
    designs: dict[str, Path] = {}
    for sheet, data in sheets.items():
        stem = "table1" if sheet.endswith("1") else "table2"
        count_path = raw_dir / f"nakao2023_{stem}_counts.tsv"
        write_count_table(data, count_path)
        counts[sheet] = count_path

    summaries = []
    commands: dict[str, dict[str, list[str]]] = {}
    for spec in ANALYSES:
        source = sheets[spec["sheet"]]
        count_path = counts[spec["sheet"]]
        count_columns = [column for column in source.columns if column not in {"sgRNA", "Gene"}]
        design = raw_dir / f"{spec['id']}.design.tsv"
        write_design_matrix(spec, count_columns, design)
        designs[spec["id"]] = design
        drugz_path, _, drugz_command = run_drugz(root, count_path, spec, raw_dir)
        mageck_path, mageck_command = run_mageck(root, count_path, design, spec, raw_dir)
        summaries.append(browser_summary(spec, source, drugz_path, mageck_path))
        commands[spec["id"]] = {"drugz": drugz_command, "mageck": mageck_command}

    output = pd.concat(summaries, ignore_index=True)
    result_path = output_dir / "gene_results.tsv"
    output.to_csv(result_path, sep="\t", index=False, na_rep="")
    manifest = {
        "study": "Nakao2023",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_workbook": str(WORKBOOK),
        "source_workbook_sha256": sha256(workbook),
        "source_sheets": EXPECTED_SHEETS,
        "browser_primary_score": "DrugZ normZ; direction-aware p/FDR/rank",
        "drugz": {"version": "1.1.0.2 build 116", "normalization": 10000000, "pseudocount": 5, "half_window_size": 500},
        "mageck": {"version": mageck_version, "normalization": "median"},
        "commands": commands,
        "result_rows": int(output.shape[0]),
        "result_sha256": sha256(result_path),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {result_path} ({output.shape[0]} rows)")


if __name__ == "__main__":
    main()
