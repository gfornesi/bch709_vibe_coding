#!/usr/bin/env python3
"""
Main objective: Analyze GC content distribution of the specified file and produce both a table and histogram as output. 

[Environment]: Write python code that runs in the [bch709_vibe_coding] conda environment. Micromamba is mapped to conda. 
view the file "copilot-instructions.md" for mroe information about installed packages and the environment. 
**Input specification:**
- File: mrna.fa.gz 
- Structure: This is a fasta file. First item in header is accession, fifth item in header is length (len)
- Additional inputs: sacCer3.fa.gz (for reference)

**Output 1 — Table:**
- Filename: please save this table in the results folder and name the file: "mrna_metrics.tsv"
- Columns: Accession, length, gc_content. Sort the values by gc_content descending 
- Decimal places: use 4 decimal places
- Sorting: Sort the values in this table by gc_content descending

**Output 2 — Plot:**
- Filename: save this plot in the results folder and name the file "gc_content_distribution.png
- Size: 1600x900 px, dpi: 200
- Colors: Please use blue for the color of the bars in the histogram and red for the mean vertical dashed lines and orange for the median vertical dashed line.
- Axes/labels/legend: This is a histogram with density curve overlays. Display the mean and median values.   Label the x-axis "GC Content" and have numbers ranging from 0-1 on the x-axis.
"""

import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

# File paths
mrna_path = Path("data/mrna.fa.gz")
output_table = Path("results/mrna_metrics.tsv")
output_plot = Path("results/gc_content_distribution.png")

# Ensure results directory exists
Path("results").mkdir(exist_ok=True)


def parse_fasta_gz(filepath):
    """
    Parse a gzipped FASTA file and extract accession, length, and sequence.
    Accession: first item in header
    Length: fifth item in header
    """
    records = []
    with gzip.open(filepath, 'rt') as f:
        current_header = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process previous record if exists
                if current_header:
                    header_parts = current_header.split()
                    accession = header_parts[0][1:]  # Remove '>'
                    length = int(header_parts[4]) if len(header_parts) > 4 else len(''.join(current_seq))
                    sequence = ''.join(current_seq)
                    records.append({
                        'header': current_header,
                        'accession': accession,
                        'length': length,
                        'sequence': sequence
                    })
                
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Process last record
        if current_header:
            header_parts = current_header.split()
            accession = header_parts[0][1:]  # Remove '>'
            length = int(header_parts[4]) if len(header_parts) > 4 else len(''.join(current_seq))
            sequence = ''.join(current_seq)
            records.append({
                'header': current_header,
                'accession': accession,
                'length': length,
                'sequence': sequence
            })
    
    return records


def calculate_gc_content(sequence):
    """Calculate GC content as a fraction (0-1)."""
    if len(sequence) == 0:
        return 0.0
    
    # Count G and C (case-insensitive)
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    gc_content = gc_count / len(sequence)
    
    return gc_content


def main():
    print(f"Reading FASTA file: {mrna_path}")
    records = parse_fasta_gz(mrna_path)
    print(f"Read {len(records)} sequences")
    
    # Calculate GC content for each sequence
    gc_data = []
    for record in records:
        gc_content = calculate_gc_content(record['sequence'])
        gc_data.append({
            'Accession': record['accession'],
            'length': record['length'],
            'gc_content': gc_content
        })
    
    # Create DataFrame
    df = pd.DataFrame(gc_data)
    
    # Sort by gc_content descending
    df = df.sort_values('gc_content', ascending=False).reset_index(drop=True)
    
    # Round gc_content to 4 decimal places
    df['gc_content'] = df['gc_content'].round(4)
    
    # Save table
    print(f"Saving table to: {output_table}")
    df.to_csv(output_table, sep='\t', index=False)
    
    print(f"\nTable summary:")
    print(f"Total sequences: {len(df)}")
    print(f"Mean GC content: {df['gc_content'].mean():.4f}")
    print(f"Median GC content: {df['gc_content'].median():.4f}")
    print(f"Min GC content: {df['gc_content'].min():.4f}")
    print(f"Max GC content: {df['gc_content'].max():.4f}")
    
    # Create histogram with density curve overlay
    print(f"\nCreating histogram: {output_plot}")
    
    fig, ax = plt.subplots(figsize=(16, 9), dpi=200)
    
    gc_values = df['gc_content'].values
    mean_gc = gc_values.mean()
    median_gc = np.median(gc_values)
    
    # Histogram with density
    counts, bins, patches = ax.hist(gc_values, bins=50, density=True, 
                                      color='blue', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add density curve
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(gc_values)
    x_range = np.linspace(0, 1, 1000)
    ax.plot(x_range, kde(x_range), 'k-', linewidth=2, label='Density')
    
    # Add mean and median lines
    ax.axvline(mean_gc, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_gc:.4f}')
    ax.axvline(median_gc, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_gc:.4f}')
    
    # Set labels and limits
    ax.set_xlabel('GC Content', fontsize=12, fontweight='bold')
    ax.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax.set_title('GC Content Distribution', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=200, bbox_inches='tight')
    print(f"Histogram saved successfully!")
    
    plt.close()


if __name__ == '__main__':
    main()
