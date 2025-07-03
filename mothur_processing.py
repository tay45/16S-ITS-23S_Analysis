#!/usr/bin/env python3
"""
mothur_processing_2_fixed.py (length‑filtered + group‑aware)

Pipeline steps
--------------
1. summary.seqs
2. screen.seqs   (filters < min_len or > max_len, default 1 000–3 000 bp)
3. classify.seqs (uses filtered *.good.fasta + *.good.groups)
4. remove.lineage (optional)

Change log
~~~~~~~~~~
2025‑07‑02 • FIX: correctly pass the group file produced by screen.seqs ("*.good.groups").
2025‑07‑03 • FIX: use the correct parameter name (`taxon=`) in remove.lineage and allow a
               plain‑text exclude file (one taxon per line) or a literal list.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

# ───────────────────────── helpers ────────────────────────────────

def setup_logging(log_file: str, level: str = "INFO") -> None:
    Path(log_file).expanduser().parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )
    logging.debug("Current PATH: %s", os.environ.get("PATH"))


def run_mothur(mothur_path: str, cmd: str, label: str) -> None:
    """Execute a mothur one‑liner and exit on failure."""
    logging.info("Running mothur %s …", label)
    logging.debug("Command: %s", cmd)
    try:
        res = subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if res.stdout:
            logging.debug(res.stdout.strip())
        if res.stderr:
            logging.warning(res.stderr.strip())
        logging.info("%s finished OK.", label)
    except subprocess.CalledProcessError as e:
        logging.error("%s failed:\n%s", label, e.stderr)
        sys.exit(1)

# ───────────────────────── mothur steps ───────────────────────────

def summary_seqs(mothur: str, fasta: str, procs: int) -> None:
    cmd = f'"{mothur}" "#summary.seqs(fasta={fasta}, processors={procs})"'
    run_mothur(mothur, cmd, "summary.seqs")


def screen_seqs(
    mothur: str,
    fasta: str,
    group: str,
    min_len: int,
    max_len: int,
    procs: int,
) -> tuple[str, str]:
    cmd = (
        f'"{mothur}" "#screen.seqs('
        f'fasta={fasta}, group={group}, minlength={min_len}, '
        f'maxlength={max_len}, maxambig=0, processors={procs})"'
    )
    run_mothur(mothur, cmd, "screen.seqs")

    good_fasta = Path(fasta).with_suffix("").as_posix() + ".good.fasta"
    good_group = Path(group).with_suffix("").as_posix() + ".good.groups"

    for pth in (good_fasta, good_group):
        if not Path(pth).is_file():
            logging.error("screen.seqs did not create %s", pth)
            sys.exit(1)

    return good_fasta, good_group


def classify_seqs(
    mothur: str,
    fasta: str,
    group: str,
    reference: str,
    taxonomy: str,
    method: str,
    numwanted: int,
    search: str,
    procs: int,
) -> None:
    cmd = (
        f'"{mothur}" "#classify.seqs('
        f'fasta={fasta}, group={group}, reference={reference}, '
        f'taxonomy={taxonomy}, method={method}, numwanted={numwanted}, '
        f'search={search}, processors={procs})"'
    )
    run_mothur(mothur, cmd, "classify.seqs")


def build_taxon_arg(lineage_exclude: str) -> str:
    """Parse a lineage‑exclude file or literal list into a mothur taxon string."""
    p = Path(lineage_exclude)
    if p.is_file():
        taxa = [ln.strip() for ln in p.read_text().splitlines() if ln.strip() and not ln.startswith('#')]
        if not taxa:
            logging.error("%s is empty — cannot build taxon list for remove.lineage", lineage_exclude)
            sys.exit(1)
        taxon_arg = "-".join(taxa)
    else:
        taxon_arg = lineage_exclude  # assume user passed literal list (e.g., "Chloroplast-Mitochondria")
    logging.debug("remove.lineage taxon argument: %s", taxon_arg)
    return taxon_arg


def remove_lineage(
    mothur: str,
    fasta: str,
    taxonomy: str,
    lineage_exclude: str,
) -> None:
    taxon_arg = build_taxon_arg(lineage_exclude)
    cmd = (
        f'"{mothur}" "#remove.lineage('
        f'fasta={fasta}, taxonomy={taxonomy}, taxon={taxon_arg})"'
    )
    run_mothur(mothur, cmd, "remove.lineage")

# ───────────────────────── path resolver ──────────────────────────

def resolve_mothur(user_path: Optional[str]) -> str:
    DEFAULT = "/opt/mothur/1.48/bin/mothur"
    mothur = user_path or shutil.which("mothur") or (DEFAULT if Path(DEFAULT).is_file() else None)
    if not mothur or not Path(mothur).is_file():
        logging.error("Cannot find mothur; pass --mothur-path or load module.")
        sys.exit(1)
    if not os.access(mothur, os.X_OK):
        logging.error("mothur at %s is not executable.", mothur)
        sys.exit(1)
    return mothur

# ─────────────────────────── main ────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser("Mothur processing pipeline")
    p.add_argument("--combined-fasta", required=True)
    p.add_argument("--combined-group", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--reference-fasta", required=True)
    p.add_argument("--taxonomy-file", required=True)
    p.add_argument("--method", default="knn")
    p.add_argument("--numwanted", type=int, default=1)
    p.add_argument("--search", default="blastplus")
    p.add_argument("--processors", type=int, default=8)
    p.add_argument("--remove-lineage", action="store_true")
    p.add_argument("--lineage-exclude", default="Chloroplast-Mitochondria")
    p.add_argument("--min-length", type=int, default=1000)
    p.add_argument("--max-length", type=int, default=3000)
    p.add_argument("--log-file", required=True)
    p.add_argument("--mothur-path")
    args = p.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    setup_logging(args.log_file, "DEBUG")

    mothur = resolve_mothur(args.mothur_path)
    logging.info("Using mothur: %s", mothur)

    summary_seqs(mothur, args.combined_fasta, args.processors)

    good_fasta, good_group = screen_seqs(
        mothur,
        fasta=args.combined_fasta,
        group=args.combined_group,
        min_len=args.min_length,
        max_len=args.max_length,
        procs=args.processors,
    )

    classify_seqs(
        mothur,
        fasta=good_fasta,
        group=good_group,
        reference=args.reference_fasta,
        taxonomy=args.taxonomy_file,
        method=args.method,
        numwanted=args.numwanted,
        search=args.search,
        procs=args.processors,
    )

    if args.remove_lineage:
        base_fa = Path(good_fasta).stem
        base_tax = Path(args.taxonomy_file).stem
        tax_out = Path(args.output_dir) / f"{base_fa}.{base_tax}.{args.method}.taxonomy"
        if not tax_out.is_file():
            logging.error("Taxonomy file not found: %s", tax_out)
            sys.exit(1)
        remove_lineage(
            mothur,
            fasta=good_fasta,
            taxonomy=str(tax_out),
            lineage_exclude=args.lineage_exclude,
        )

    logging.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
