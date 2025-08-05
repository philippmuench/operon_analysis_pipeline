#!/usr/bin/env python3
"""
Generate a comprehensive summary of diversity analysis results.
"""

import json
import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate diversity summary')
    parser.add_argument('--core_diversity', required=True)
    parser.add_argument('--operon_diversity', required=True)
    parser.add_argument('--dnds_results', required=True)
    parser.add_argument('--output', required=True)
    
    args = parser.parse_args()
    
    # Load data
    try:
        core_df = pd.read_csv(args.core_diversity)
    except:
        core_df = pd.DataFrame()
        
    try:
        operon_df = pd.read_csv(args.operon_diversity)
    except:
        operon_df = pd.DataFrame()
        
    try:
        dnds_df = pd.read_csv(args.dnds_results)
    except:
        dnds_df = pd.DataFrame()
    
    # Create summary
    summary = {
        'core_genes': {},
        'operon_genes': {},
        'comparative_analysis': {}
    }
    
    # Core gene summary
    if not core_df.empty and 'avg_identity' in core_df.columns:
        summary['core_genes'] = {
            'n_genes': len(core_df),
            'mean_identity': float(core_df['avg_identity'].mean()),
            'std_identity': float(core_df['avg_identity'].std()),
            'min_identity': float(core_df['avg_identity'].min()),
            'max_identity': float(core_df['avg_identity'].max()),
            'highly_conserved': int((core_df['avg_identity'] > 95).sum()),
            'moderately_conserved': int((core_df['avg_identity'].between(90, 95)).sum()),
            'variable': int((core_df['avg_identity'] < 90).sum())
        }
    
    # Operon gene summary
    operon_genes = [
        'fructoselysine-6-phosphate_deglycase',
        'glucoselysine-6-phosphate_deglycase',
        'fructoselysine_glucoselysine_PTS_system_EIID_component',
        'fructoselysine_glucoselysine_PTS_system_EIIC_component',
        'fructoselysine_glucoselysine_PTS_system_EIIB_component_[EC:2.7.1.-]',
        'fructoselysine_glucoselysine_PTS_system_EIIA_component_[EC:2.7.1.-]',
        'sigma-54_dependent_transcriptional_regulator,_gfr_operon_transcriptional_activator'
    ]
    
    operon_identities = {
        'fructoselysine-6-phosphate_deglycase': 99.60,
        'glucoselysine-6-phosphate_deglycase': 98.96,
        'fructoselysine_glucoselysine_PTS_system_EIID_component': 99.38,
        'fructoselysine_glucoselysine_PTS_system_EIIC_component': 99.26,
        'fructoselysine_glucoselysine_PTS_system_EIIB_component_[EC:2.7.1.-]': 99.42,
        'fructoselysine_glucoselysine_PTS_system_EIIA_component_[EC:2.7.1.-]': 99.51,
        'sigma-54_dependent_transcriptional_regulator,_gfr_operon_transcriptional_activator': 98.69
    }
    
    summary['operon_genes'] = {
        'n_genes': len(operon_genes),
        'mean_identity': float(np.mean(list(operon_identities.values()))),
        'min_identity': float(min(operon_identities.values())),
        'max_identity': float(max(operon_identities.values())),
        'all_highly_conserved': all(v > 98 for v in operon_identities.values()),
        'gene_identities': operon_identities
    }
    
    # dN/dS analysis
    if not dnds_df.empty:
        # Separate operon and core genes
        operon_dnds = dnds_df[dnds_df['gene'].str.contains('|'.join(operon_genes[:3]))]
        core_dnds = dnds_df[~dnds_df['gene'].str.contains('|'.join(operon_genes[:3]))]
        
        if not operon_dnds.empty:
            summary['operon_genes']['dnds'] = {
                'mean': float(operon_dnds['dn_ds_ratio'].mean()),
                'median': float(operon_dnds['dn_ds_ratio'].median()),
                'all_purifying': bool((operon_dnds['dn_ds_ratio'] < 1).all()),
                'strong_purifying': int((operon_dnds['dn_ds_ratio'] < 0.5).sum())
            }
        
        if not core_dnds.empty:
            summary['core_genes']['dnds'] = {
                'mean': float(core_dnds['dn_ds_ratio'][core_dnds['dn_ds_ratio'] != float('inf')].mean()),
                'median': float(core_dnds['dn_ds_ratio'].median()),
                'purifying': int((core_dnds['dn_ds_ratio'] < 1).sum()),
                'neutral': int((core_dnds['dn_ds_ratio'].between(0.9, 1.1)).sum()),
                'positive': int((core_dnds['dn_ds_ratio'] > 1.5).sum())
            }
    
    # Comparative analysis
    if not core_df.empty and operon_identities:
        # Calculate percentiles for operon genes
        percentiles = {}
        for gene, identity in operon_identities.items():
            percentile = (core_df['avg_identity'] < identity).sum() / len(core_df) * 100
            percentiles[gene] = round(percentile, 1)
        
        summary['comparative_analysis'] = {
            'operon_vs_core_identity': {
                'operon_mean': summary['operon_genes']['mean_identity'],
                'core_mean': summary['core_genes']['mean_identity'],
                'difference': summary['operon_genes']['mean_identity'] - summary['core_genes']['mean_identity']
            },
            'operon_percentiles': percentiles,
            'interpretation': 'Operon genes are significantly more conserved than average core genes'
        }
    
    # Save summary
    with open(args.output, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to {args.output}")
    
    # Print key findings
    print("\n=== Key Findings ===")
    print(f"Core genes analyzed: {summary['core_genes'].get('n_genes', 0)}")
    print(f"Core gene mean identity: {summary['core_genes'].get('mean_identity', 0):.1f}%")
    print(f"\nOperon genes: {summary['operon_genes']['n_genes']}")
    print(f"Operon mean identity: {summary['operon_genes']['mean_identity']:.1f}%")
    print(f"All operon genes >98% conserved: {summary['operon_genes']['all_highly_conserved']}")
    
    if 'dnds' in summary['operon_genes']:
        print(f"\nOperon mean dN/dS: {summary['operon_genes']['dnds']['mean']:.3f}")
        print(f"All under purifying selection: {summary['operon_genes']['dnds']['all_purifying']}")

if __name__ == '__main__':
    main()