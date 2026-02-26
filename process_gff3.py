#!/usr/bin/env python3
"""
Process GFF3 file to count genes, exons, tRNA, and snoRNA per chromosome.
Output feature counts and densities per Mb.
"""

import gzip
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

# File paths
gff3_path = Path("data/saccharomyces_cerevisiae.gff.gz")
chrom_sizes_path = Path("data/chrom.sizes")
output_tsv = Path("results/chr_feature_counts.tsv")
dropped_seqids_file = Path("results/dropped_seqids.txt")

# Ensure results directory exists
Path("results").mkdir(exist_ok=True)

# Read chromosome sizes
chrom_sizes = {}
with open(chrom_sizes_path, 'r') as f:
    for line in f:
        chrom, length = line.strip().split('\t')
        chrom_sizes[chrom] = int(length)

# Initialize counters for each chromosome
counters = defaultdict(lambda: {
    'n_gene': 0,
    'n_exon_unique': set(),  # Store unique (start, end, strand) tuples
    'n_tRNA': 0,
    'n_snoRNA': 0
})

# Initialize all chromosomes with 0 counts
for chrom in chrom_sizes.keys():
    counters[chrom] = {
        'n_gene': 0,
        'n_exon_unique': set(),
        'n_tRNA': 0,
        'n_snoRNA': 0
    }

# Track dropped seqids
dropped_seqids = set()
excluded_feature_lines = 0

# Read GFF3 file (gzipped)
with gzip.open(gff3_path, 'rt') as f:
    for line in f:
        # Skip header lines
        if line.startswith('#'):
            continue
        
        # Parse GFF3 line
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        seqid = fields[0]
        gff_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        
        # Check if seqid is in chrom.sizes
        if seqid not in chrom_sizes:
            dropped_seqids.add(seqid)
            excluded_feature_lines += 1
            continue
        
        # Count features by type
        if gff_type == "gene":
            counters[seqid]['n_gene'] += 1
        elif gff_type in ("exon", "noncoding_exon", "CDS"):
            # Store unique (start, end, strand) tuples for exons/CDS
            # This counts both protein-coding (CDS) and non-coding (noncoding_exon) exon-like features
            counters[seqid]['n_exon_unique'].add((start, end, strand))
        elif gff_type == "tRNA":
            counters[seqid]['n_tRNA'] += 1
        elif gff_type == "snoRNA":
            counters[seqid]['n_snoRNA'] += 1

# Convert sets to counts
for chrom in counters:
    counters[chrom]['n_exon_unique'] = len(counters[chrom]['n_exon_unique'])

# Create output dataframe
data = []
for chrom in chrom_sizes.keys():
    c = counters[chrom]
    data.append({
        'chrom': chrom,
        'chrom_length_bp': chrom_sizes[chrom],
        'n_gene': c['n_gene'],
        'n_exon_unique': c['n_exon_unique'],
        'n_tRNA': c['n_tRNA'],
        'n_snoRNA': c['n_snoRNA']
    })

df = pd.DataFrame(data)

# Calculate densities per Mb
df['gene_per_Mb'] = df['n_gene'] / (df['chrom_length_bp'] / 1e6)
df['exon_unique_per_Mb'] = df['n_exon_unique'] / (df['chrom_length_bp'] / 1e6)
df['tRNA_per_Mb'] = df['n_tRNA'] / (df['chrom_length_bp'] / 1e6)
df['snoRNA_per_Mb'] = df['n_snoRNA'] / (df['chrom_length_bp'] / 1e6)

# Round to 4 decimal places
df['gene_per_Mb'] = df['gene_per_Mb'].round(4)
df['exon_unique_per_Mb'] = df['exon_unique_per_Mb'].round(4)
df['tRNA_per_Mb'] = df['tRNA_per_Mb'].round(4)
df['snoRNA_per_Mb'] = df['snoRNA_per_Mb'].round(4)

# Sort by gene_per_Mb descending
df = df.sort_values('gene_per_Mb', ascending=False).reset_index(drop=True)

# Save to TSV
df.to_csv(output_tsv, sep='\t', index=False)

# Save dropped seqids
dropped_seqids_sorted = sorted(dropped_seqids)
with open(dropped_seqids_file, 'w') as f:
    for seqid in dropped_seqids_sorted:
        f.write(f"{seqid}\n")

# Print summary to console
print(f"Number of dropped seqids: {len(dropped_seqids_sorted)}")
print(f"Number of excluded feature lines: {excluded_feature_lines}")
print("\nTop 5 rows of results:")
print(df.head().to_string(index=False))

# Create visualizations
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Gene counts per chromosome
ax1 = axes[0, 0]
df_plot = df.sort_values('n_gene', ascending=True)
ax1.barh(df_plot['chrom'], df_plot['n_gene'], color='steelblue')
ax1.set_xlabel('Number of Genes')
ax1.set_title('Gene Counts per Chromosome')
ax1.grid(axis='x', alpha=0.3)

# Plot 2: Gene density per Mb
ax2 = axes[0, 1]
df_plot2 = df.sort_values('gene_per_Mb', ascending=True)
ax2.barh(df_plot2['chrom'], df_plot2['gene_per_Mb'], color='coral')
ax2.set_xlabel('Gene Density (per Mb)')
ax2.set_title('Gene Density per Chromosome')
ax2.grid(axis='x', alpha=0.3)

# Plot 3: Feature counts across all chromosomes
ax3 = axes[1, 0]
df_sorted = df.sort_values('chrom')
ax3.plot(df_sorted['chrom'], df_sorted['n_gene'], marker='o', label='Genes', linewidth=2)
ax3.plot(df_sorted['chrom'], df_sorted['n_exon_unique'], marker='s', label='Exons', linewidth=2)
ax3.plot(df_sorted['chrom'], df_sorted['n_tRNA'], marker='^', label='tRNA', linewidth=2)
ax3.plot(df_sorted['chrom'], df_sorted['n_snoRNA'], marker='d', label='snoRNA', linewidth=2)
ax3.set_ylabel('Feature Count')
ax3.set_title('Feature Counts Across Chromosomes')
ax3.legend()
ax3.grid(alpha=0.3)
ax3.tick_params(axis='x', rotation=45)

# Plot 4: Stacked bar chart of density per Mb
ax4 = axes[1, 1]
df_sorted = df.sort_values('chrom')
width = 0.6
ax4.bar(df_sorted['chrom'], df_sorted['gene_per_Mb'], width, label='Gene', color='steelblue')
ax4.bar(df_sorted['chrom'], df_sorted['exon_unique_per_Mb'], width, 
        bottom=df_sorted['gene_per_Mb'], label='Exon', color='coral')
ax4.bar(df_sorted['chrom'], df_sorted['tRNA_per_Mb'], width, 
        bottom=df_sorted['gene_per_Mb'] + df_sorted['exon_unique_per_Mb'], label='tRNA', color='lightgreen')
ax4.bar(df_sorted['chrom'], df_sorted['snoRNA_per_Mb'], width, 
        bottom=df_sorted['gene_per_Mb'] + df_sorted['exon_unique_per_Mb'] + df_sorted['tRNA_per_Mb'], 
        label='snoRNA', color='gold')
ax4.set_ylabel('Density (per Mb)')
ax4.set_title('Feature Density Composition per Chromosome')
ax4.legend(loc='upper right')
ax4.tick_params(axis='x', rotation=45)
ax4.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('results/feature_counts_visualization.png', dpi=300, bbox_inches='tight')
print("\nVisualization saved to results/feature_counts_visualization.png")
