#!/usr/bin/env python3
"""Record the four prepared caches without reading or rebuilding raw inputs."""
import csv
import hashlib
from pathlib import Path

CACHES = (
    "gene_scores.rds", "depmap_crispr.rds", "external_context.rds",
    "additional_screen_candidates.rds",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "data/metadata/cache_manifest.csv"
    rows = []
    for name in CACHES:
        path = root / "data/processed" / name
        rows.append({"path": path.relative_to(root).as_posix(),
                     "bytes": path.stat().st_size, "sha256": sha256(path)})
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("path", "bytes", "sha256"), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Recorded {len(rows)} caches in {output.relative_to(root)}")


if __name__ == "__main__":
    main()
