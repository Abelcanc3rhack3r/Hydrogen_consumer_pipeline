#!/usr/bin/env python3
"""
Create a QIIME 2 paired-end FASTQ manifest from a reads directory.

Usage:
    python ./make_manifest.py ./reads ./manifest_mmdb1.csv

Expected read-name examples:
    sample1_R1.fastq.gz / sample1_R2.fastq.gz
    sample1_1.fastq.gz  / sample1_2.fastq.gz
    sample1.forward.fastq.gz / sample1.reverse.fastq.gz

The output manifest uses QIIME 2's paired-end manifest columns:
    sample-id,absolute-filepath,direction

For portability, the filepath values are written relative to the current
working directory instead of as machine-specific absolute paths.
"""

import csv
import re
import sys
from pathlib import Path

FASTQ_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)

FORWARD_PATTERNS = [
    re.compile(r"(.+?)(?:[_\.-]R?1)(?:[_\.-].*)?$", re.IGNORECASE),
    re.compile(r"(.+?)(?:[_\.-]forward)(?:[_\.-].*)?$", re.IGNORECASE),
]

REVERSE_PATTERNS = [
    re.compile(r"(.+?)(?:[_\.-]R?2)(?:[_\.-].*)?$", re.IGNORECASE),
    re.compile(r"(.+?)(?:[_\.-]reverse)(?:[_\.-].*)?$", re.IGNORECASE),
]


def strip_fastq_suffix(filename: str) -> str:
    lower = filename.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def infer_sample_and_direction(path: Path):
    stem = strip_fastq_suffix(path.name)

    for pattern in FORWARD_PATTERNS:
        match = pattern.match(stem)
        if match:
            return clean_sample_id(match.group(1)), "forward"

    for pattern in REVERSE_PATTERNS:
        match = pattern.match(stem)
        if match:
            return clean_sample_id(match.group(1)), "reverse"

    return None, None


def clean_sample_id(sample_id: str) -> str:
    sample_id = sample_id.strip()
    sample_id = re.sub(r"[_\.-]+$", "", sample_id)
    return sample_id


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python ./make_manifest.py ./reads ./manifest_mmdb1.csv", file=sys.stderr)
        return 1

    reads_dir = Path(sys.argv[1])
    output_csv = Path(sys.argv[2])

    if not reads_dir.exists() or not reads_dir.is_dir():
        print(f"ERROR: reads directory not found: {reads_dir}", file=sys.stderr)
        return 1

    fastq_files = sorted(
        path for path in reads_dir.iterdir()
        if path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)
    )

    if not fastq_files:
        print(f"ERROR: no FASTQ files found in {reads_dir}", file=sys.stderr)
        return 1

    samples = {}
    skipped = []

    for path in fastq_files:
        sample_id, direction = infer_sample_and_direction(path)
        if sample_id is None:
            skipped.append(path.name)
            continue
        samples.setdefault(sample_id, {})[direction] = path

    complete_samples = {
        sample_id: directions
        for sample_id, directions in samples.items()
        if "forward" in directions and "reverse" in directions
    }

    incomplete_samples = {
        sample_id: directions
        for sample_id, directions in samples.items()
        if sample_id not in complete_samples
    }

    if not complete_samples:
        print("ERROR: no complete paired-end samples were found.", file=sys.stderr)
        print("Expected names like sample_R1.fastq.gz and sample_R2.fastq.gz.", file=sys.stderr)
        return 1

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample-id", "absolute-filepath", "direction"])

        for sample_id in sorted(complete_samples):
            for direction in ("forward", "reverse"):
                path = complete_samples[sample_id][direction]
                relative_path = Path(path).as_posix()
                writer.writerow([sample_id, relative_path, direction])

    print(f"Wrote manifest for {len(complete_samples)} paired-end samples: {output_csv}")

    if incomplete_samples:
        print("WARNING: incomplete samples were skipped:", file=sys.stderr)
        for sample_id, directions in sorted(incomplete_samples.items()):
            found = ", ".join(sorted(directions))
            print(f"  {sample_id}: found {found}", file=sys.stderr)

    if skipped:
        print("WARNING: files with unrecognized read direction were skipped:", file=sys.stderr)
        for filename in skipped:
            print(f"  {filename}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
