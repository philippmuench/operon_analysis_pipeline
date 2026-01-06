#!/usr/bin/env python3
"""
Filter BLAST output files to remove 7 excluded assemblies.
Run from 03_blast_search directory.
"""

import os
import pandas as pd
from datetime import datetime

# Assemblies to exclude
EXCLUDE_PATTERNS = [
    'ENT_CA1913AA_AS', 'ENT_CA1914AA_AS', 'ENT_CA1915AA_AS',
    'ENT_CA1916AA_AS', 'ENT_CA1917AA_AS', 'ENT_CA1918AA_AS', 'ENT_CA1919AA_AS'
]

def should_exclude(value):
    """Check if a value contains any excluded assembly pattern."""
    if pd.isna(value):
        return False
    value_str = str(value)
    return any(pattern in value_str for pattern in EXCLUDE_PATTERNS)


def filter_blast_csv(filepath):
    """Filter a BLAST CSV file to remove excluded genomes."""
    print(f"\nFiltering: {filepath}")
    df = pd.read_csv(filepath)
    original_count = len(df)

    # Try different column names for genome ID
    genome_cols = ['genome_id', 'genome', 'sseqid', 'qseqid']
    filtered = False

    for col in genome_cols:
        if col in df.columns:
            mask = ~df[col].apply(should_exclude)
            df_filtered = df[mask]
            if len(df_filtered) < original_count:
                df = df_filtered
                filtered = True
                print(f"  Filtered on column '{col}': {original_count} -> {len(df)} rows")
                break

    if not filtered:
        print(f"  No matching column found, checking all columns...")
        for col in df.columns:
            if df[col].dtype == 'object':
                mask = ~df[col].apply(should_exclude)
                if mask.sum() < len(df):
                    df = df[mask]
                    print(f"  Filtered on column '{col}': {original_count} -> {len(df)} rows")
                    break

    df.to_csv(filepath, index=False)
    print(f"  Saved: {len(df)} rows")
    return original_count - len(df)


def main():
    print("=" * 60)
    print("Filtering BLAST Output Files to 8,596 Genomes")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # BLAST CSV files to filter (relative to 03_blast_search)
    blast_files = [
        'output/all_blast_hits.csv',
        'output/all_blast_hits_detailed.csv',
        'output/ptsA_filtered_hits.csv',
        'output/ptsB_filtered_hits.csv',
        'output/ptsC_filtered_hits.csv',
        'output/ptsD_filtered_hits.csv',
    ]

    total_removed = 0
    for filepath in blast_files:
        if os.path.exists(filepath):
            removed = filter_blast_csv(filepath)
            total_removed += removed
        else:
            print(f"\nSkipping (not found): {filepath}")

    print(f"\n" + "=" * 60)
    print(f"Total rows removed: {total_removed}")
    print("DONE!")
    print("=" * 60)


if __name__ == '__main__':
    main()
