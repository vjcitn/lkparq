#!/usr/bin/env python3
"""
GTF to Parquet Converter - Full Attribute Preservation

Converts GENCODE/Ensembl GTF files to Parquet format, preserving ALL attributes
unlike the lossy TxDb schema.

Output structure:
    output_dir/
    ├── genes.parquet          # One row per gene
    ├── transcripts.parquet    # One row per transcript
    ├── exons.parquet          # One row per exon (per transcript)
    ├── cds.parquet            # CDS regions
    ├── features.parquet       # All features (start/stop codons, UTRs, etc.)
    └── metadata.parquet       # GTF header metadata

Usage:
    python gtf_to_parquet.py input.gtf output_directory/
"""

import argparse
import re
import gzip
from pathlib import Path
from collections import defaultdict
from typing import Iterator

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# =============================================================================
# GTF Parsing
# =============================================================================

def parse_attributes(attr_string: str) -> dict:
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    # Match: key "value" or key value
    pattern = r'(\S+)\s+"([^"]*)"'
    
    for match in re.finditer(pattern, attr_string):
        key, value = match.groups()
        # Handle multiple tags by collecting into list
        if key in attrs:
            if isinstance(attrs[key], list):
                attrs[key].append(value)
            else:
                attrs[key] = [attrs[key], value]
        else:
            attrs[key] = value
    
    return attrs


def parse_gtf_line(line: str) -> dict | None:
    """Parse a single GTF line into a dictionary."""
    if line.startswith('#'):
        return None
    
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 9:
        return None
    
    chrom, source, feature, start, end, score, strand, frame, attributes = fields
    
    record = {
        'chrom': chrom,
        'source': source,
        'feature': feature,
        'start': int(start),
        'end': int(end),
        'score': score if score != '.' else None,
        'strand': strand,
        'frame': frame if frame != '.' else None,
    }
    
    # Parse attributes
    attrs = parse_attributes(attributes)
    record.update(attrs)
    
    return record


def parse_gtf_header(filepath: str) -> dict:
    """Extract metadata from GTF header comments."""
    metadata = {}
    
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    with open_func(filepath, mode) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('##'):
                # Parse ##key: value format
                match = re.match(r'##(\S+):\s*(.+)', line.strip())
                if match:
                    metadata[match.group(1)] = match.group(2)
    
    # Extract genome build from description if present
    # Handles patterns like "human genome (GRCh38)" or "mouse genome (GRCm39)"
    if 'description' in metadata:
        genome_match = re.search(r'\(([A-Za-z0-9_.]+)\)', metadata['description'])
        if genome_match:
            metadata['genome'] = genome_match.group(1)
    
    # Also check for explicit genome-build field (some GTFs have this)
    if 'genome-build' in metadata and 'genome' not in metadata:
        metadata['genome'] = metadata['genome-build']
    
    return metadata


def iter_gtf(filepath: str) -> Iterator[dict]:
    """Iterate over GTF records."""
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    with open_func(filepath, mode) as f:
        for line in f:
            record = parse_gtf_line(line)
            if record is not None:
                yield record


# =============================================================================
# Data Collection
# =============================================================================

def collect_features(gtf_path: str) -> dict[str, list]:
    """
    Collect all features from GTF, organized by type.
    
    Returns dict with keys: genes, transcripts, exons, cds, other_features
    """
    genes = {}  # gene_id -> gene record
    transcripts = {}  # transcript_id -> transcript record
    exons = []
    cds_regions = []
    other_features = []  # start_codon, stop_codon, UTR, etc.
    
    print("  Parsing GTF file...")
    record_count = 0
    
    for record in iter_gtf(gtf_path):
        record_count += 1
        if record_count % 500000 == 0:
            print(f"    Processed {record_count:,} records...")
        
        feature_type = record['feature']
        
        if feature_type == 'gene':
            gene_id = record.get('gene_id', '')
            genes[gene_id] = {
                'gene_id': gene_id,
                'gene_id_stripped': gene_id.split('.')[0],
                'gene_name': record.get('gene_name'),
                'gene_type': record.get('gene_type'),
                'chrom': record['chrom'],
                'start': record['start'],
                'end': record['end'],
                'strand': record['strand'],
                'source': record['source'],
                'level': _safe_int(record.get('level')),
                'tags': _normalize_tags(record.get('tag')),
                'havana_gene': record.get('havana_gene'),
            }
        
        elif feature_type == 'transcript':
            tx_id = record.get('transcript_id', '')
            transcripts[tx_id] = {
                'transcript_id': tx_id,
                'transcript_id_stripped': tx_id.split('.')[0],
                'transcript_name': record.get('transcript_name'),
                'transcript_type': record.get('transcript_type'),
                'gene_id': record.get('gene_id', ''),
                'gene_id_stripped': record.get('gene_id', '').split('.')[0],
                'gene_name': record.get('gene_name'),
                'chrom': record['chrom'],
                'start': record['start'],
                'end': record['end'],
                'strand': record['strand'],
                'source': record['source'],
                'level': _safe_int(record.get('level')),
                'tags': _normalize_tags(record.get('tag')),
                'transcript_support_level': record.get('transcript_support_level'),
                'havana_transcript': record.get('havana_transcript'),
                'ccdsid': record.get('ccdsid'),
                'protein_id': record.get('protein_id'),
            }
        
        elif feature_type == 'exon':
            exons.append({
                'exon_id': record.get('exon_id', ''),
                'exon_id_stripped': record.get('exon_id', '').split('.')[0],
                'exon_number': _safe_int(record.get('exon_number')),
                'transcript_id': record.get('transcript_id', ''),
                'transcript_id_stripped': record.get('transcript_id', '').split('.')[0],
                'gene_id': record.get('gene_id', ''),
                'gene_id_stripped': record.get('gene_id', '').split('.')[0],
                'chrom': record['chrom'],
                'start': record['start'],
                'end': record['end'],
                'strand': record['strand'],
                'source': record['source'],
                'level': _safe_int(record.get('level')),
                'tags': _normalize_tags(record.get('tag')),
            })
        
        elif feature_type == 'CDS':
            cds_regions.append({
                'cds_id': f"{record.get('transcript_id', '')}:CDS:{record['start']}-{record['end']}",
                'transcript_id': record.get('transcript_id', ''),
                'transcript_id_stripped': record.get('transcript_id', '').split('.')[0],
                'gene_id': record.get('gene_id', ''),
                'gene_id_stripped': record.get('gene_id', '').split('.')[0],
                'protein_id': record.get('protein_id'),
                'exon_number': _safe_int(record.get('exon_number')),
                'chrom': record['chrom'],
                'start': record['start'],
                'end': record['end'],
                'strand': record['strand'],
                'frame': _safe_int(record.get('frame')),
                'source': record['source'],
                'level': _safe_int(record.get('level')),
                'tags': _normalize_tags(record.get('tag')),
            })
        
        elif feature_type in ('start_codon', 'stop_codon', 'UTR', 
                               'five_prime_utr', 'three_prime_utr', 'Selenocysteine'):
            other_features.append({
                'feature_type': feature_type,
                'transcript_id': record.get('transcript_id', ''),
                'transcript_id_stripped': record.get('transcript_id', '').split('.')[0],
                'gene_id': record.get('gene_id', ''),
                'gene_id_stripped': record.get('gene_id', '').split('.')[0],
                'chrom': record['chrom'],
                'start': record['start'],
                'end': record['end'],
                'strand': record['strand'],
                'frame': _safe_int(record.get('frame')),
                'source': record['source'],
                'level': _safe_int(record.get('level')),
                'tags': _normalize_tags(record.get('tag')),
            })
    
    print(f"  Total records: {record_count:,}")
    print(f"    Genes: {len(genes):,}")
    print(f"    Transcripts: {len(transcripts):,}")
    print(f"    Exons: {len(exons):,}")
    print(f"    CDS: {len(cds_regions):,}")
    print(f"    Other features: {len(other_features):,}")
    
    return {
        'genes': list(genes.values()),
        'transcripts': list(transcripts.values()),
        'exons': exons,
        'cds': cds_regions,
        'other_features': other_features,
    }


def _safe_int(value) -> int | None:
    """Safely convert to int, returning None on failure."""
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _normalize_tags(tag_value) -> str | None:
    """Normalize tag(s) to semicolon-separated string."""
    if tag_value is None:
        return None
    if isinstance(tag_value, list):
        return ';'.join(tag_value)
    return tag_value


# =============================================================================
# Schema Definitions
# =============================================================================

GENE_SCHEMA = pa.schema([
    ('gene_id', pa.string()),
    ('gene_id_stripped', pa.string()),
    ('gene_name', pa.string()),
    ('gene_type', pa.string()),
    ('chrom', pa.string()),
    ('start', pa.int64()),
    ('end', pa.int64()),
    ('strand', pa.string()),
    ('source', pa.string()),
    ('level', pa.int32()),
    ('tags', pa.string()),
    ('havana_gene', pa.string()),
])

TRANSCRIPT_SCHEMA = pa.schema([
    ('transcript_id', pa.string()),
    ('transcript_id_stripped', pa.string()),
    ('transcript_name', pa.string()),
    ('transcript_type', pa.string()),
    ('gene_id', pa.string()),
    ('gene_id_stripped', pa.string()),
    ('gene_name', pa.string()),
    ('chrom', pa.string()),
    ('start', pa.int64()),
    ('end', pa.int64()),
    ('strand', pa.string()),
    ('source', pa.string()),
    ('level', pa.int32()),
    ('tags', pa.string()),
    ('transcript_support_level', pa.string()),
    ('havana_transcript', pa.string()),
    ('ccdsid', pa.string()),
    ('protein_id', pa.string()),
])

EXON_SCHEMA = pa.schema([
    ('exon_id', pa.string()),
    ('exon_id_stripped', pa.string()),
    ('exon_number', pa.int32()),
    ('transcript_id', pa.string()),
    ('transcript_id_stripped', pa.string()),
    ('gene_id', pa.string()),
    ('gene_id_stripped', pa.string()),
    ('chrom', pa.string()),
    ('start', pa.int64()),
    ('end', pa.int64()),
    ('strand', pa.string()),
    ('source', pa.string()),
    ('level', pa.int32()),
    ('tags', pa.string()),
])

CDS_SCHEMA = pa.schema([
    ('cds_id', pa.string()),
    ('transcript_id', pa.string()),
    ('transcript_id_stripped', pa.string()),
    ('gene_id', pa.string()),
    ('gene_id_stripped', pa.string()),
    ('protein_id', pa.string()),
    ('exon_number', pa.int32()),
    ('chrom', pa.string()),
    ('start', pa.int64()),
    ('end', pa.int64()),
    ('strand', pa.string()),
    ('frame', pa.int32()),
    ('source', pa.string()),
    ('level', pa.int32()),
    ('tags', pa.string()),
])

FEATURE_SCHEMA = pa.schema([
    ('feature_type', pa.string()),
    ('transcript_id', pa.string()),
    ('transcript_id_stripped', pa.string()),
    ('gene_id', pa.string()),
    ('gene_id_stripped', pa.string()),
    ('chrom', pa.string()),
    ('start', pa.int64()),
    ('end', pa.int64()),
    ('strand', pa.string()),
    ('frame', pa.int32()),
    ('source', pa.string()),
    ('level', pa.int32()),
    ('tags', pa.string()),
])


# =============================================================================
# Parquet Export
# =============================================================================

def write_parquet(df: pd.DataFrame, path: Path, schema: pa.Schema) -> None:
    """Write DataFrame to Parquet with schema."""
    table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
    pq.write_table(
        table,
        path,
        compression='snappy',
        write_statistics=True,
    )


def export_to_parquet(
    features: dict[str, list],
    metadata: dict,
    output_dir: Path,
    partition_genes_by_chrom: bool = False,
) -> None:
    """Export all features to Parquet files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Export genes
    if features['genes']:
        print("  Writing genes.parquet...")
        df = pd.DataFrame(features['genes'])
        
        if partition_genes_by_chrom:
            table = pa.Table.from_pandas(df, schema=GENE_SCHEMA, preserve_index=False)
            pq.write_to_dataset(
                table,
                root_path=str(output_dir / 'genes'),
                partition_cols=['chrom'],
                compression='snappy',
            )
        else:
            df = df.sort_values(['chrom', 'start'])
            write_parquet(df, output_dir / 'genes.parquet', GENE_SCHEMA)
    
    # Export transcripts
    if features['transcripts']:
        print("  Writing transcripts.parquet...")
        df = pd.DataFrame(features['transcripts']).sort_values(['chrom', 'start'])
        write_parquet(df, output_dir / 'transcripts.parquet', TRANSCRIPT_SCHEMA)
    
    # Export exons
    if features['exons']:
        print("  Writing exons.parquet...")
        df = pd.DataFrame(features['exons']).sort_values(['chrom', 'start'])
        write_parquet(df, output_dir / 'exons.parquet', EXON_SCHEMA)
    
    # Export CDS
    if features['cds']:
        print("  Writing cds.parquet...")
        df = pd.DataFrame(features['cds']).sort_values(['chrom', 'start'])
        write_parquet(df, output_dir / 'cds.parquet', CDS_SCHEMA)
    
    # Export other features
    if features['other_features']:
        print("  Writing features.parquet...")
        df = pd.DataFrame(features['other_features']).sort_values(['chrom', 'start'])
        write_parquet(df, output_dir / 'features.parquet', FEATURE_SCHEMA)
    
    # Export metadata
    print("  Writing metadata.parquet...")
    meta_records = [{'name': k, 'value': v} for k, v in metadata.items()]
    meta_df = pd.DataFrame(meta_records)
    meta_schema = pa.schema([('name', pa.string()), ('value', pa.string())])
    write_parquet(meta_df, output_dir / 'metadata.parquet', meta_schema)


# =============================================================================
# Main
# =============================================================================

def convert_gtf_to_parquet(
    input_gtf: str,
    output_dir: str,
    partition_by_chrom: bool = False,
) -> None:
    """
    Convert GTF file to Parquet format with full attribute preservation.
    
    Args:
        input_gtf: Path to input GTF file (can be gzipped)
        output_dir: Output directory for Parquet files
        partition_by_chrom: Whether to partition genes by chromosome
    """
    input_path = Path(input_gtf)
    output_path = Path(output_dir)
    
    if not input_path.exists():
        raise FileNotFoundError(f"GTF file not found: {input_path}")
    
    print(f"Parsing {input_path}...")
    
    # Parse header metadata
    metadata = parse_gtf_header(str(input_path))
    print(f"  Metadata: {metadata}")
    
    # Collect all features
    features = collect_features(str(input_path))
    
    # Export to Parquet
    print(f"\nWriting Parquet files to {output_path}...")
    export_to_parquet(features, metadata, output_path, partition_by_chrom)
    
    print("\nConversion complete!")
    print(f"Output files in: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GTF to Parquet with full attribute preservation"
    )
    parser.add_argument(
        "input_gtf",
        help="Input GTF file (can be .gz)"
    )
    parser.add_argument(
        "output_dir",
        help="Output directory for Parquet files"
    )
    parser.add_argument(
        "--partition-by-chrom",
        action="store_true",
        help="Partition gene output by chromosome"
    )
    
    args = parser.parse_args()
    
    convert_gtf_to_parquet(
        args.input_gtf,
        args.output_dir,
        args.partition_by_chrom,
    )


if __name__ == "__main__":
    main()
