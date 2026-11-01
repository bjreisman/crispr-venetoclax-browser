#!/usr/bin/env python3
"""Copy the explicit public file set to a new checkout; leave local Git intact.

Usage: python3 scripts/export_publication.py /absolute/path/to/new-directory
Review and commit that directory with Git before pushing. Raw data and existing
Git history are deliberately excluded. Dependencies: Python standard library.
"""
import argparse
from pathlib import Path
import shutil

from update_cache_manifest import sha256


def public_files(root):
    fixed = [
        ".gitignore", ".gitattributes", ".rscignore", "README.md", "LICENSE", "CITATION.cff",
        "app.R", "setup.R", "renv.lock", "data/gene_sets_curated.csv",
        "scripts/requirements-nakao2023.txt", "scripts/vendor/drugz.py",
        "scripts/vendor/DRUGZ_LICENSE",
    ]
    files = {root / name for name in fixed}
    for directory, patterns in (
        ("R", ("*.R",)), ("scripts", ("*.R", "*.py")),
        ("docs", ("*.md",)), ("data/metadata", ("*.csv", "*.txt", "*.json")),
        (".github/workflows", ("*.yml", "*.yaml")),
    ):
        for pattern in patterns:
            files.update((root / directory).glob(pattern))
    return sorted(files)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = args.destination.resolve()
    if destination.exists():
        parser.error("Destination already exists; choose a new empty path.")
    if destination == root or root in destination.parents:
        parser.error("Destination must be outside the research project.")
    files = public_files(root)
    for path in files:
        if path.is_symlink() or not path.is_file():
            parser.error(f"Required regular file missing or symlinked: {path}")
        if path.stat().st_size >= 100 * 1024 * 1024:
            parser.error(f"File exceeds GitHub's ordinary Git limit: {path}")
    destination.mkdir(parents=True)
    total = 0
    for path in files:
        relative = path.relative_to(root)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        if sha256(path) != sha256(target):
            raise RuntimeError(f"Copy checksum mismatch: {relative}")
        total += target.stat().st_size
        print(relative.as_posix())
    print(f"Exported {len(files)} files ({total / 1024**2:.1f} MiB) to {destination}")


if __name__ == "__main__":
    main()
