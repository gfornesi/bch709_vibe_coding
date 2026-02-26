# GFF3 Feature Counter

This repository contains a Python script to process a GFF3 file of *Saccharomyces cerevisiae* features, count various feature types per chromosome, and generate summary tables and plots.

## Files

- `data/saccharomyces_cerevisiae.gff.gz` – input GFF3 (compressed) with genomic annotations.
- `data/chrom.sizes` – chromosome lengths (TSV, two columns: `chrom` and `length_bp`).
- `results/chr_feature_counts.tsv` – output table with feature counts and densities.
- `results/dropped_seqids.txt` – list of seqids from GFF that were not in `chrom.sizes`.
- `results/feature_counts_visualization.png` – generated summary plots.
- `process_gff3.py` – main Python processing script.

## Requirements

The script is designed to run within the `bch709` conda/micromamba environment. The following Python packages must be installed:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `biopython` (not directly used in current script but available)
- `tqdm` (not currently used; included for future enhancements)

You can install these via micromamba/conda:

```bash
micromamba install -n bch709 pandas numpy matplotlib seaborn biopython tqdm -y
```
## Usage

Run the script from the repository root with the appropriate environment active:

```bash
# activate environment (example with micromamba)
micromamba run -n bch709 python3 process_gff3.py
```

The script will:

1. Read chromosome sizes from `data/chrom.sizes`.
2. Parse the GFF3 file, counting:
   - genes (`type == gene`)
   - unique exon/CDS intervals (`type in exon, noncoding_exon, CDS`), deduped by `(start,end,strand)`
   - tRNA (`type == tRNA`)
   - snoRNA (`type == snoRNA`)
3. Exclude any features whose seqid is not listed in `chrom.sizes`; these seqids are logged in `results/dropped_seqids.txt` and the count of excluded lines printed.
4. Build a summary table with densities per megabase, sorted by gene density.
5. Save results to `results/chr_feature_counts.tsv` and produce a four-panel visualization (`results/feature_counts_visualization.png`).
6. Print a brief report to the console including the number of dropped seqids, excluded lines, and the top 5 rows of the result table.

## Output Format

The TSV has the following columns:

- `chrom`
- `chrom_length_bp`
- `n_gene`
- `n_exon_unique`
- `n_tRNA`
- `n_snoRNA`
- `gene_per_Mb`
- `exon_unique_per_Mb`
- `tRNA_per_Mb`
- `snoRNA_per_Mb`

Densities are rounded to four decimal places. All chromosomes from `chrom.sizes` appear, with zeroes for missing features.

## Notes

- The script currently treats `chrmt` (mitochondrial) as dropped since it isn't listed in `chrom.sizes`.
- Visualizations are saved automatically; feel free to modify or extend plotting code in `process_gff3.py`.

## Manual

1. Ensure data files are in place under `data/`.
2. Activate or run within the `bch709` environment.
3. Execute `python process_gff3.py`.
4. Examine the outputs in `results/`.
5. To add new chromosome sizes or annotation files, update `data/chrom.sizes` or change the script accordingly.

---

This project was created by the GitHub Copilot assistant on behalf of the user.

