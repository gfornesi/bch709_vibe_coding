# Project Instructions for Copilot

This repository contains a Python script to process a GFF3 file of *Saccharomyces cerevisiae* features, count various feature types per chromosome, and generate summary tables and plots.

## Files

- `data/saccharomyces_cerevisiae.gff.gz` – input GFF3 (compressed) with genomic annotations.
- `data/chrom.sizes` – chromosome lengths (TSV, two columns: `chrom` and `length_bp`).
- `results/chr_feature_counts.tsv` – output table with feature counts and densities.
- `results/dropped_seqids.txt` – list of seqids from GFF that were not in `chrom.sizes`.
- `results/feature_counts_visualization.png` – generated summary plots.
- `process_gff3.py` – main Python processing script.

## Environment

The code runs inside a `micromamba`/`conda` environment named `bch709`. Packages installed:

- pandas
- numpy
- matplotlib
- seaborn
- biopython
- tqdm

The shell is Linux-based. Use micromamba commands such as:

```bash
micromamba env list        # to view environments
micromamba run -n bch709 python3 process_gff3.py
micromamba install -n bch709 pandas numpy matplotlib seaborn biopython tqdm -y
```

If conda activation is needed, initialize the shell with `eval "$(micromamba shell hook --shell bash)"` and then `micromamba activate bch709`.

## Usage

Run the script from the project root:

```bash
micromamba run -n bch709 python3 process_gff3.py
```

Outputs: counts TSV and visualization PNG in `results/`, plus a report printed to console.

## Notes

- Only chromosomes listed in `data/chrom.sizes` are processed; others are logged to `results/dropped_seqids.txt`.
- Input annotations are filtered by seqid, and counts are generated per chromosome.
- Exon counting is based on unique (start,end,strand) intervals across types `exon`, `noncoding_exon`, and `CDS`.

This file aims to help GitHub Copilot understand repository purpose, environment setup, and usage instructions.