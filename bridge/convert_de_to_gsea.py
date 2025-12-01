#!/usr/bin/env python3
"""
Bridge Script: Convert DE Analysis Results to GSEA Pipeline Format

This script converts the output from RNA-Seq_DE_GO_analysis pipeline
to the input format required by RNA-Seq_GO_GSEA_analysis pipeline.

Input: final_de_results.csv from DE analysis
Output: Excel file compatible with GSEA pipeline

Author: Pipeline Integration Team
Date: 2025-12-01
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
import yaml
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_de_results(csv_path):
    """
    Load DE analysis results from CSV file.
    
    Args:
        csv_path: Path to final_de_results.csv
        
    Returns:
        pandas.DataFrame
    """
    logger.info(f"Loading DE results from: {csv_path}")
    
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"  Loaded {len(df)} genes")
        logger.info(f"  Columns: {', '.join(df.columns)}")
        return df
    except Exception as e:
        logger.error(f"Failed to load CSV: {e}")
        raise


def standardize_column_names(df, source_pipeline="deseq2"):
    """
    Standardize column names to match GSEA pipeline expectations.
    
    The GSEA pipeline expects:
    - Gene ID column (any name containing 'gene', 'id', or 'symbol')
    - log2FoldChange or logFC
    - pvalue
    - padj or adj.P.Val
    - baseMean or AveExpr
    
    Args:
        df: Input DataFrame
        source_pipeline: Source pipeline type ("deseq2", "edger", "limma")
        
    Returns:
        DataFrame with standardized column names
    """
    logger.info("Standardizing column names...")
    
    df_copy = df.copy()
    
    # Column mapping based on source pipeline
    # Most columns are already standardized from DE pipeline
    column_mapping = {
        # DESeq2 format (already good)
        'log2FoldChange': 'log2FoldChange',
        'pvalue': 'pvalue',
        'padj': 'padj',
        'baseMean': 'baseMean',
        
        # edgeR format
        'logFC': 'log2FoldChange',
        'PValue': 'pvalue',
        'FDR': 'padj',
        'logCPM': 'baseMean',
        
        # limma format
        'adj.P.Val': 'padj',
        'P.Value': 'pvalue',
        'AveExpr': 'baseMean',
        't': 'statistic'
    }
    
    # Rename columns if they exist
    rename_dict = {}
    for old_name, new_name in column_mapping.items():
        if old_name in df_copy.columns and old_name != new_name:
            rename_dict[old_name] = new_name
    
    if rename_dict:
        df_copy = df_copy.rename(columns=rename_dict)
        logger.info(f"  Renamed columns: {rename_dict}")
    
    return df_copy


def add_required_columns(df):
    """
    Add or compute required columns for GSEA pipeline.
    
    Args:
        df: Input DataFrame
        
    Returns:
        DataFrame with all required columns
    """
    logger.info("Adding required columns...")
    
    # Ensure we have a Gene column (first column is usually gene ID)
    if 'Gene' not in df.columns and 'GeneID' not in df.columns:
        # Use the first column as Gene ID
        first_col = df.columns[0]
        if first_col not in ['log2FoldChange', 'pvalue', 'padj']:
            df = df.rename(columns={first_col: 'GeneID'})
            logger.info(f"  Renamed '{first_col}' to 'GeneID'")
    
    # Add ranking metric for GSEA (if not present)
    if 'rank_metric' not in df.columns:
        if 'log2FoldChange' in df.columns and 'pvalue' in df.columns:
            # Signed ranking: sign(log2FC) * -log10(pvalue)
            import numpy as np
            df['rank_metric'] = np.sign(df['log2FoldChange']) * -np.log10(df['pvalue'].replace(0, 1e-300))
            logger.info("  Added rank_metric column")
    
    # Add regulation direction
    if 'regulation' not in df.columns and 'log2FoldChange' in df.columns:
        df['regulation'] = df['log2FoldChange'].apply(
            lambda x: 'up' if x > 0 else 'down' if x < 0 else 'no change'
        )
        logger.info("  Added regulation column")
    
    return df


def filter_and_sort(df, sort_by='padj'):
    """
    Filter invalid entries and sort by significance.
    
    Args:
        df: Input DataFrame
        sort_by: Column to sort by (default: 'padj')
        
    Returns:
        Filtered and sorted DataFrame
    """
    logger.info("Filtering and sorting...")
    
    # Remove rows with NA in critical columns
    critical_cols = ['log2FoldChange', 'pvalue', 'padj']
    before_count = len(df)
    
    for col in critical_cols:
        if col in df.columns:
            df = df.dropna(subset=[col])
    
    after_count = len(df)
    if before_count != after_count:
        logger.info(f"  Removed {before_count - after_count} rows with NA values")
    
    # Sort by padj (most significant first)
    if sort_by in df.columns:
        df = df.sort_values(sort_by)
        logger.info(f"  Sorted by {sort_by}")
    
    return df


def save_to_excel(df, output_path, comparison_name):
    """
    Save DataFrame to Excel file with proper formatting.
    
    Args:
        df: Input DataFrame
        output_path: Output Excel file path
        comparison_name: Name of the comparison (for sheet name)
    """
    logger.info(f"Saving to Excel: {output_path}")
    
    try:
        # Create Excel writer
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            # Write main results
            df.to_excel(writer, sheet_name='DE_Results', index=False)
            
            # Write metadata sheet
            metadata = pd.DataFrame({
                'Parameter': [
                    'Comparison',
                    'Total Genes',
                    'Significant Genes (padj < 0.05)',
                    'Up-regulated',
                    'Down-regulated',
                    'Source Pipeline',
                    'Conversion Date'
                ],
                'Value': [
                    comparison_name,
                    len(df),
                    len(df[df['padj'] < 0.05]) if 'padj' in df.columns else 'N/A',
                    len(df[(df['log2FoldChange'] > 0) & (df['padj'] < 0.05)]) if 'log2FoldChange' in df.columns and 'padj' in df.columns else 'N/A',
                    len(df[(df['log2FoldChange'] < 0) & (df['padj'] < 0.05)]) if 'log2FoldChange' in df.columns and 'padj' in df.columns else 'N/A',
                    'RNA-Seq_DE_GO_analysis',
                    pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
                ]
            })
            metadata.to_excel(writer, sheet_name='Metadata', index=False)
            
            # Write significant genes only (for quick reference)
            if 'padj' in df.columns:
                sig_genes = df[df['padj'] < 0.05]
                if len(sig_genes) > 0:
                    sig_genes.to_excel(writer, sheet_name='Significant_Only', index=False)
                    logger.info(f"  Added sheet 'Significant_Only' with {len(sig_genes)} genes")
        
        logger.info(f"✓ Successfully saved to {output_path}")
        
    except Exception as e:
        logger.error(f"Failed to save Excel: {e}")
        raise


def convert_de_to_gsea_format(
    input_csv,
    output_excel,
    comparison_name=None,
    source_pipeline="deseq2"
):
    """
    Main conversion function.
    
    Args:
        input_csv: Path to final_de_results.csv
        output_excel: Path to output Excel file
        comparison_name: Name of comparison (extracted from path if None)
        source_pipeline: Source pipeline type
    """
    # Extract comparison name from path if not provided
    if comparison_name is None:
        # Try to extract from path: .../pairwise/H2O2_vs_Control/final_de_results.csv
        parts = Path(input_csv).parts
        for i, part in enumerate(parts):
            if part == 'pairwise' and i + 1 < len(parts):
                comparison_name = parts[i + 1]
                break
        if comparison_name is None:
            comparison_name = "comparison"
    
    logger.info("="*60)
    logger.info(f"Converting DE results to GSEA format")
    logger.info(f"  Comparison: {comparison_name}")
    logger.info("="*60)
    
    # Load data
    df = load_de_results(input_csv)
    
    # Standardize columns
    df = standardize_column_names(df, source_pipeline)
    
    # Add required columns
    df = add_required_columns(df)
    
    # Filter and sort
    df = filter_and_sort(df)
    
    # Save to Excel
    save_to_excel(df, output_excel, comparison_name)
    
    logger.info("="*60)
    logger.info("✓ Conversion completed successfully!")
    logger.info("="*60)


def batch_convert(de_output_dir, gsea_input_dir, pattern="pairwise/*/final_de_results.csv"):
    """
    Batch convert all DE results in a directory.
    
    Args:
        de_output_dir: Root directory of DE analysis output (e.g., output/H2O2_Neuron)
        gsea_input_dir: Target directory for GSEA pipeline input
        pattern: Glob pattern to find DE result files (relative to de_output_dir)
    """
    de_output_path = Path(de_output_dir)
    gsea_input_path = Path(gsea_input_dir)
    
    # Create output directory if it doesn't exist
    gsea_input_path.mkdir(parents=True, exist_ok=True)
    
    # Find all DE result files
    result_files = list(de_output_path.glob(pattern))
    
    if not result_files:
        logger.warning(f"No files found matching pattern: {pattern}")
        return
    
    logger.info(f"Found {len(result_files)} DE result files to convert")
    logger.info("")
    
    for i, csv_file in enumerate(result_files, 1):
        # Extract comparison name from path
        parts = csv_file.parts
        comparison_name = "comparison"
        for j, part in enumerate(parts):
            if part == 'pairwise' and j + 1 < len(parts):
                comparison_name = parts[j + 1]
                break
        
        # Create output filename
        output_file = gsea_input_path / f"{comparison_name}_DE_results.xlsx"
        
        logger.info(f"[{i}/{len(result_files)}] Processing: {comparison_name}")
        
        try:
            convert_de_to_gsea_format(
                input_csv=str(csv_file),
                output_excel=str(output_file),
                comparison_name=comparison_name
            )
        except Exception as e:
            logger.error(f"  ✗ Failed: {e}")
            continue
        
        logger.info("")
    
    logger.info("="*60)
    logger.info(f"✓ Batch conversion completed!")
    logger.info(f"  Converted files saved to: {gsea_input_path}")
    logger.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Convert DE analysis results to GSEA pipeline format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert single file
  python convert_de_to_gsea.py \\
    --input output/H2O2_Neuron/pairwise/H2O2_vs_Control/final_de_results.csv \\
    --output ../RNA-Seq_GO_GSEA_analysis/data/H2O2_vs_Control.xlsx
  
  # Batch convert all comparisons
  python convert_de_to_gsea.py \\
    --batch \\
    --de-output-dir output/H2O2_Neuron \\
    --gsea-input-dir ../RNA-Seq_GO_GSEA_analysis/data/from_de_pipeline
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        help='Input CSV file (final_de_results.csv)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Output Excel file path'
    )
    
    parser.add_argument(
        '--comparison-name', '-n',
        type=str,
        help='Name of the comparison (auto-detected from path if not provided)'
    )
    
    parser.add_argument(
        '--source-pipeline', '-s',
        type=str,
        choices=['deseq2', 'edger', 'limma'],
        default='deseq2',
        help='Source DE pipeline (default: deseq2)'
    )
    
    parser.add_argument(
        '--batch', '-b',
        action='store_true',
        help='Batch mode: convert all DE results in directory'
    )
    
    parser.add_argument(
        '--de-output-dir',
        type=str,
        help='DE analysis output directory (for batch mode)'
    )
    
    parser.add_argument(
        '--gsea-input-dir',
        type=str,
        help='GSEA pipeline input directory (for batch mode)'
    )
    
    args = parser.parse_args()
    
    # Batch mode
    if args.batch:
        if not args.de_output_dir or not args.gsea_input_dir:
            parser.error("Batch mode requires --de-output-dir and --gsea-input-dir")
        
        batch_convert(args.de_output_dir, args.gsea_input_dir)
    
    # Single file mode
    else:
        if not args.input or not args.output:
            parser.error("Single file mode requires --input and --output")
        
        convert_de_to_gsea_format(
            input_csv=args.input,
            output_excel=args.output,
            comparison_name=args.comparison_name,
            source_pipeline=args.source_pipeline
        )


if __name__ == '__main__':
    main()
