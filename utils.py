"""
Utility functions for Variant Visualizer
Contains badge creation, filtering, and data processing functions
UPDATED WITH GENE NAME MAPPING AND OMIM LINKS - HANDLES MULTIPLE GENES
UPDATED WITH LARGER FONTS FOR BETTER READABILITY
SIMPLIFIED: REMOVED GREEN GENE FILTERING FROM APPLY_FILTERS
UPDATED WITH NEW GNOMAD AND CGEN FREQUENCY HANDLING
"""

import dash_bootstrap_components as dbc
import polars as pl
import pandas as pd
from dash import html
import os
from config import *

# Global variable to store gene mapping
_gene_mapping = None

def load_gene_mapping():
    """Load HGNC ID to gene name mapping from file - OPTIMIZED VERSION"""
    global _gene_mapping
    
    if _gene_mapping is not None:
        return _gene_mapping
    
    gene_mapping_file = os.path.join(DATA_DIR, 'symbol.tsv')
    
    try:
        if os.path.exists(gene_mapping_file):
            import time
            start_time = time.time()
            
            # OPTIMIZATION: Use Polars for faster loading if file is large
            try:
                import polars as pl
                mapping_df = pl.read_csv(gene_mapping_file, separator='\t')
                # Convert to dictionary efficiently
                _gene_mapping = dict(zip(
                    mapping_df['HGNC_ID'].cast(pl.Utf8).to_list(),
                    mapping_df['HGNC_symbol'].to_list()
                ))
                logger.info(f"Loaded {len(_gene_mapping)} gene mappings with Polars in {time.time() - start_time:.2f}s")
            except Exception:
                # Fallback to pandas if Polars fails
                mapping_df = pd.read_csv(gene_mapping_file, sep='\t')
                _gene_mapping = dict(zip(mapping_df['HGNC_ID'].astype(str), mapping_df['HGNC_symbol']))
                logger.info(f"Loaded {len(_gene_mapping)} gene mappings with Pandas in {time.time() - start_time:.2f}s")
        else:
            logger.warning(f"Gene mapping file not found: {gene_mapping_file}")
            _gene_mapping = {}
    except Exception as e:
        logger.error(f"Error loading gene mapping: {e}")
        _gene_mapping = {}
    
    return _gene_mapping

def get_gene_name_from_id(gene_id):
    """Convert HGNC ID to gene name - IMPROVED VERSION"""
    if pd.isna(gene_id) or gene_id in [None, '', 'UNKNOWN']:
        return gene_id if gene_id else 'UNKNOWN'
    
    gene_mapping = load_gene_mapping()
    gene_id_str = str(gene_id)
    
    # First check if it's already a gene name (not an ID)
    if gene_mapping and gene_id_str not in gene_mapping:
        # Check if the value is already a gene name
        for gid, gname in gene_mapping.items():
            if gname.upper() == gene_id_str.upper():
                return gene_id_str  # It's already a gene name
    
    # Return gene name if found in mapping, otherwise return original ID
    return gene_mapping.get(gene_id_str, gene_id_str)

def is_gene_name_or_id(value):
    """Determine if a value is a gene name or gene ID"""
    if pd.isna(value) or value in [None, '', 'UNKNOWN']:
        return 'unknown'
    
    gene_mapping = load_gene_mapping()
    value_str = str(value)
    
    # Check if it's a gene ID (found as key in mapping)
    if value_str in gene_mapping:
        return 'gene_id'
    
    # Check if it's a gene name (found as value in mapping)
    for gene_id, gene_name in gene_mapping.items():
        if gene_name.upper() == value_str.upper():
            return 'gene_name'
    
    # Check if it looks like an ID (numeric or starts with HGNC)
    if value_str.isdigit() or value_str.startswith('HGNC:'):
        return 'potential_gene_id'
    
    # Assume it's a gene name if it's alphabetic
    if value_str.isalpha():
        return 'potential_gene_name'
    
    return 'unknown'

def find_genes_by_search_term(search_term):
    """Find all genes (IDs and names) that match a search term"""
    if not search_term or not search_term.strip():
        return [], []
    
    search_lower = search_term.strip().lower()
    gene_mapping = load_gene_mapping()
    
    matching_ids = []
    matching_names = []
    
    if not gene_mapping:
        logger.warning("No gene mapping available for search")
        return matching_ids, matching_names
    
    # Search in both gene IDs and gene names
    for gene_id, gene_name in gene_mapping.items():
        # Search in gene name (primary search)
        if search_lower in gene_name.lower():
            matching_ids.append(gene_id)
            matching_names.append(gene_name)
        # Also search in gene ID (for partial ID matches)
        elif search_lower in gene_id.lower():
            matching_ids.append(gene_id)
            matching_names.append(gene_name)
    
    return matching_ids, matching_names

def debug_search_data(df, search_term="BRCA1"):
    """Debug function to understand your data structure"""
    print(f"\n=== DEBUG: Search Data Analysis ===")
    
    if isinstance(df, pl.DataFrame):
        print(f"DataFrame type: Polars, Shape: {df.shape}")
        if len(df) > 0:
            print(f"Columns: {df.columns}")
            sample_row = df.row(0, named=True)
        else:
            print("DataFrame is empty")
            return
    else:
        print(f"DataFrame type: Pandas, Shape: {df.shape}")
        if len(df) > 0:
            print(f"Columns: {list(df.columns)}")
            sample_row = df.iloc[0].to_dict()
        else:
            print("DataFrame is empty")
            return
    
    print(f"\nSample row data:")
    for key in ['CHROM', 'POS', 'gene', 'SAMPLE', 'consequence']:
        if key in sample_row:
            print(f"  {key}: {sample_row[key]} (type: {type(sample_row[key])})")
    
    # Check gene mapping
    gene_mapping = load_gene_mapping()
    print(f"\nGene mapping: {len(gene_mapping)} entries loaded")
    
    if gene_mapping:
        # Test search
        matching_ids, matching_names = find_genes_by_search_term(search_term)
        print(f"Search for '{search_term}': {len(matching_ids)} matches")
        if matching_names:
            print(f"Found genes: {matching_names[:5]}")
    
    # Check unique values in gene column
    if isinstance(df, pl.DataFrame):
        unique_genes = df.select('gene').unique().limit(10).to_pandas()['gene'].tolist()
    else:
        unique_genes = df['gene'].unique()[:10].tolist()
    
    print(f"\nSample gene values in data: {unique_genes}")

def apply_filters_optimized(df, search_term=None, genotype_filter=None, chromosome_filter=None, 
                          active_filters=None, selected_samples=None):
    """OPTIMIZED apply filters for large datasets - MUCH FASTER"""
    
    if len(df) == 0:
        logger.info("Input DataFrame is empty")
        return df
    
    # Convert to Polars if it's pandas for faster operations
    if isinstance(df, pd.DataFrame):
        df = pl.from_pandas(df)
    
    logger.info(f"Starting optimized filtering with {len(df)} variants")
    
    # Sample filter (most selective first)
    if selected_samples:
        df = df.filter(pl.col('SAMPLE').is_in(selected_samples))
        logger.info(f"After sample filter: {len(df)} variants")
    
    # Chromosome filter
    if chromosome_filter and chromosome_filter != "all":
        df = df.filter(pl.col('CHROM') == chromosome_filter)
        logger.info(f"After chromosome filter: {len(df)} variants")
    
    # Genotype filter
    if genotype_filter and genotype_filter != "all":
        if genotype_filter == "het":
            df = df.filter(pl.col('GT').is_in(['0/1', '1/0', '0|1', '1|0']))
        elif genotype_filter == "hom_alt":
            df = df.filter(pl.col('GT').is_in(['1/1', '1|1']))
        elif genotype_filter == "hom_ref":
            df = df.filter(pl.col('GT').is_in(['0/0', '0|0']))
        logger.info(f"After genotype filter: {len(df)} variants")
    
    # OPTIMIZED SEARCH FILTER - Much faster for large datasets
    if search_term and search_term.strip():
        search_lower = search_term.strip().lower()
        logger.info(f"Applying OPTIMIZED search filter: '{search_term}'")
        
        search_conditions = []
        
        # OPTIMIZATION: Use vectorized string operations
        try:
            # Sample search - vectorized
            search_conditions.append(pl.col('SAMPLE').str.to_lowercase().str.contains(search_lower))
        except Exception as e:
            logger.error(f"Sample search error: {e}")
        
        # Consequence search - vectorized
        try:
            search_conditions.append(pl.col('consequence').str.to_lowercase().str.contains(search_lower))
        except Exception as e:
            logger.error(f"Consequence search error: {e}")
        
        # Chromosome search - optimized
        try:
            chrom_search = search_lower.replace('chr', '').replace('chromosome', '').strip()
            if chrom_search:
                search_conditions.extend([
                    pl.col('CHROM').str.to_lowercase() == chrom_search,
                    pl.col('CHROM').str.to_lowercase().str.contains(chrom_search)
                ])
        except Exception as e:
            logger.error(f"Chromosome search error: {e}")
        
        # Position search - optimized
        try:
            if search_term.isdigit():
                pos_value = int(search_term)
                search_conditions.append(pl.col('POS') == pos_value)
            
            # Handle chr:pos format
            if ':' in search_term:
                parts = search_term.split(':')
                if len(parts) >= 2:
                    chrom_part = parts[0].replace('chr', '').strip().lower()
                    pos_part = parts[1].strip()
                    if pos_part.isdigit():
                        pos_value = int(pos_part)
                        search_conditions.append(
                            (pl.col('CHROM').str.to_lowercase() == chrom_part) & 
                            (pl.col('POS') == pos_value)
                        )
        except Exception as e:
            logger.error(f"Position search error: {e}")
        
        # OPTIMIZED Gene search - Pre-compute gene sets for faster lookup
        try:
            # Direct gene field search
            search_conditions.append(pl.col('gene').str.to_lowercase().str.contains(search_lower))
            
            # OPTIMIZATION: Pre-compute gene mapping once
            gene_mapping = load_gene_mapping()
            if gene_mapping:
                matching_gene_ids, matching_gene_names = find_genes_by_search_term(search_term)
                
                if matching_gene_ids:
                    search_conditions.append(pl.col('gene').is_in(matching_gene_ids))
                    
                if matching_gene_names:
                    search_conditions.append(pl.col('gene').is_in(matching_gene_names))
                    
        except Exception as e:
            logger.error(f"Gene search error: {e}")
        
        # Apply combined search if we have conditions - OPTIMIZED
        if search_conditions:
            # OPTIMIZATION: Combine conditions efficiently
            combined_search = search_conditions[0]
            for condition in search_conditions[1:]:
                combined_search = combined_search | condition
            
            df_before = len(df)
            df = df.filter(combined_search)
            df_after = len(df)
            
            logger.info(f"OPTIMIZED search filter result: {df_before} -> {df_after} variants")
        else:
            logger.warning(f"No valid search conditions generated for: '{search_term}'")
    
    # OPTIMIZED Preset filters - Use vectorized operations
    if active_filters:
        if active_filters.get('high_impact'):
            high_impact = ['frameshift_variant', 'stop_gained', 'stopgain', 'stop_lost']
            df = df.filter(pl.col('consequence').is_in(high_impact))
            logger.info(f"After high_impact filter: {len(df)} variants")
        
        if active_filters.get('pathogenic'):
            df = df.filter(pl.col('clinvar_sig').is_in(['Pathogenic', 'Likely pathogenic']))
            logger.info(f"After pathogenic filter: {len(df)} variants")
        
        if active_filters.get('heterozygous'):
            df = df.filter(pl.col('GT').is_in(['0/1', '1/0', '0|1', '1|0']))
            logger.info(f"After heterozygous filter: {len(df)} variants")
        
        if active_filters.get('homozygous'):
            df = df.filter(pl.col('GT').is_in(['1/1', '1|1', '0/0', '0|0']))
            logger.info(f"After homozygous filter: {len(df)} variants")
    
    logger.info(f"OPTIMIZED filtering complete: {len(df)} variants")
    return df

# Keep the original function for backward compatibility
def apply_filters(df, search_term=None, genotype_filter=None, chromosome_filter=None, 
                 active_filters=None, selected_samples=None, panel_genes=None):
    """Apply filters - SIMPLIFIED (removed green gene specific logic)"""
    return apply_filters_optimized(df, search_term, genotype_filter, chromosome_filter, 
                                 active_filters, selected_samples)

def create_gene_link(gene_id_or_name):
    """Create clickable gene name(s) with OMIM links - handles multiple genes with deduplication - LARGER FONTS"""
    
    if pd.isna(gene_id_or_name) or gene_id_or_name in [None, '', 'UNKNOWN']:
        return html.Span("UNKNOWN", style={"color": "#6c757d", "fontStyle": "italic", "fontSize": "14px"})
    
    gene_input = str(gene_id_or_name)
    
    # Handle multiple genes separated by common delimiters
    gene_separators = [' • ', ' · ', ',', ';', '|', '/', '&']
    genes = [gene_input]
    
    # Split by each separator
    for sep in gene_separators:
        if sep in gene_input:
            genes = [g.strip() for part in genes for g in part.split(sep) if g.strip()]
            break
    
    # DEDUPLICATION: Remove duplicates while preserving order
    unique_genes = []
    seen = set()
    for gene in genes:
        gene = gene.strip()
        if gene and gene not in seen:
            unique_genes.append(gene)
            seen.add(gene)
    
    # Convert each gene ID to gene name and create links
    gene_links = []
    
    for i, gene in enumerate(unique_genes):
        gene_name = get_gene_name_from_id(gene)
        
        if gene_name == 'UNKNOWN' or not gene_name:
            # If we can't map it, show the original value
            gene_element = html.Span(
                gene, 
                style={
                    "color": "#6c757d", 
                    "fontStyle": "italic",
                    "fontSize": "14px"  # INCREASED FONT SIZE
                }
            )
        else:
            # Create OMIM search URL
            omim_url = f"https://www.omim.org/search?index=entry&start=1&limit=10&search={gene_name}"
            
            gene_element = html.A(
                gene_name,
                href=omim_url,
                target="_blank",
                style={
                    "fontWeight": "bold", 
                    "color": "#0097A7", 
                    "textDecoration": "none",
                    "fontSize": "14px"  # INCREASED FONT SIZE
                },
                title=f"View {gene_name} in OMIM database",
                className="gene-link"
            )
        
        gene_links.append(gene_element)
        
        # Add separator between genes (but not after the last one)
        if i < len(unique_genes) - 1:
            gene_links.append(
                html.Span(" • ", style={
                    "color": "#6c757d", 
                    "fontSize": "14px",  # INCREASED FONT SIZE
                    "margin": "0 2px"
                })
            )
    
    # Return single element or container with multiple genes
    if len(gene_links) == 1:
        return gene_links[0]
    else:
        return html.Div(
            gene_links,
            style={
                "display": "flex",
                "flexWrap": "wrap",
                "alignItems": "center",
                "gap": "2px"
            }
        )

def is_dataframe_empty(df):
    """Check if DataFrame is empty (works for both Polars and Pandas)"""
    if isinstance(df, pl.DataFrame):
        return len(df) == 0
    else:
        return df.empty

def get_consequence_badge(consequence):
    """Return styled badge for consequence - LARGER FONTS"""
    if not consequence or pd.isna(consequence):
        consequence = 'variant'
    
    color = CONSEQUENCE_COLORS.get(consequence, 'secondary')
    return dbc.Badge(consequence, color=color, className="me-1", 
                    style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE

def get_clinvar_badge(classification):
    """Return styled badge for ClinVar classification - LARGER FONTS"""
    if not classification or pd.isna(classification):
        return dbc.Badge("No annotation", color="light", text_color="dark", 
                        style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE
    
    color = CLINVAR_COLORS.get(classification, 'secondary')
    return dbc.Badge(classification, color=color, 
                    style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE

def get_status_badge(status):
    """Return styled badge for review status - LARGER FONTS"""
    if not status or pd.isna(status):
        status = 'Pending'
    
    color = STATUS_COLORS.get(status, 'secondary')
    return dbc.Badge(status, color=color, 
                    style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE

def get_genotype_badge(genotype):
    """Return styled badge for genotype - LARGER FONTS"""
    if pd.isna(genotype) or genotype in [None, './.', '.', '']:
        return html.Span("./.", className="genotype-badge gt-missing", 
                        style={"fontSize": "0.9em", "padding": "4px 8px", "fontWeight": "600"})  # INCREASED SIZE

    gt_str = str(genotype) if genotype is not None else "./."
    gt_str = gt_str.replace('|', '/')

    base_style = {"fontSize": "0.9em", "padding": "4px 8px", "fontWeight": "600"}  # INCREASED SIZE

    if gt_str in ['0/0', '0|0']:
        return html.Span("0/0", className="genotype-badge gt-hom-ref", title="Homozygous Reference", style=base_style)
    elif gt_str in ['0/1', '1/0', '0|1', '1|0']:
        return html.Span("0/1", className="genotype-badge gt-het", title="Heterozygous", style=base_style)
    elif gt_str in ['1/1', '1|1']:
        return html.Span("1/1", className="genotype-badge gt-hom-alt", title="Homozygous Alternate", style=base_style)
    else:
        return html.Span(gt_str, className="genotype-badge gt-missing", style=base_style)

# OPTIMIZED VERSIONS FOR PERFORMANCE - WITH LARGER FONTS
def get_genotype_badge_optimized(genotype):
    """Optimized genotype badge with minimal DOM operations - LARGER FONTS"""
    if pd.isna(genotype) or genotype in [None, './.', '.', '']:
        return html.Span("./.", className="genotype-badge gt-missing")

    gt_str = str(genotype).replace('|', '/')
    
    # Use dictionary to avoid multiple if statements
    gt_classes = {
        '0/0': "genotype-badge gt-hom-ref",
        '0/1': "genotype-badge gt-het", 
        '1/0': "genotype-badge gt-het",
        '1/1': "genotype-badge gt-hom-alt"
    }
    
    return html.Span(gt_str, className=gt_classes.get(gt_str, "genotype-badge gt-missing"))

def get_consequence_badge_optimized(consequence):
    """Optimized consequence badge - LARGER FONTS"""
    if not consequence or pd.isna(consequence):
        consequence = 'variant'
    
    color = CONSEQUENCE_COLORS.get(consequence, 'secondary')
    return dbc.Badge(consequence, color=color, className="me-1", 
                    style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE

def get_clinvar_badge_optimized(classification):
    """Optimized ClinVar badge - LARGER FONTS"""
    if not classification or pd.isna(classification):
        return dbc.Badge("No annotation", color="light", text_color="dark", 
                        style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE
    
    color = CLINVAR_COLORS.get(classification, 'secondary')
    return dbc.Badge(classification, color=color, 
                    style={"fontSize": "0.8em", "padding": "0.4em 0.6em"})  # INCREASED SIZE

def format_frequency(frequency):
    """Format allele frequency for display with improved precision"""
    if pd.isna(frequency) or frequency == 0 or frequency is None:
        return "N/A"
    
    try:
        freq_val = float(frequency)
        
        if freq_val == 0:
            return "0"
        elif freq_val < 0.00001:  # < 0.001%
            return f"{freq_val:.2e}"
        elif freq_val < 0.0001:   # < 0.01%
            return f"{freq_val:.1e}"
        elif freq_val < 0.001:    # < 0.1%
            return f"{freq_val:.5f}"
        elif freq_val < 0.01:     # < 1%
            return f"{freq_val:.4f}"
        elif freq_val < 0.1:      # < 10%
            return f"{freq_val:.3f}"
        else:                     # >= 10%
            return f"{freq_val:.2f}"
            
    except (ValueError, TypeError):
        return "N/A"

def get_frequency_color_style(frequency):
    """Get color styling based on frequency value for population AF display"""
    if pd.isna(frequency) or frequency == 0 or frequency is None:
        return {"color": "#6c757d"}  # Gray for N/A
    
    try:
        freq_val = float(frequency)
        
        if freq_val < 0.001:      # < 0.1% - Very rare (red)
            return {"color": "#dc3545", "fontWeight": "bold"}
        elif freq_val < 0.01:     # 0.1% - 1% - Rare (orange)
            return {"color": "#fd7e14", "fontWeight": "bold"}
        elif freq_val < 0.05:     # 1% - 5% - Uncommon (yellow)
            return {"color": "#ffc107", "fontWeight": "bold"}
        else:                     # >= 5% - Common (green)
            return {"color": "#28a745", "fontWeight": "bold"}
            
    except (ValueError, TypeError):
        return {"color": "#6c757d"}

def format_percentage(value):
    """Format value as percentage"""
    if pd.isna(value) or value == 0:
        return "N/A"
    return f"{value:.1%}"

def format_score(score):
    """Format prediction scores for display"""
    if pd.isna(score) or score == 0:
        return "N/A"
    return f"{score:.3f}"

def get_max_gnomad_af_from_variant(variant):
    """Get the maximum gnomAD AF (excluding ASJ/FIN) from a variant record"""
    try:
        # Check if max_gnomad_af is already calculated
        if 'max_gnomad_af' in variant and pd.notna(variant['max_gnomad_af']) and variant['max_gnomad_af'] > 0:
            return variant['max_gnomad_af']
        
        # Calculate from population-specific frequencies
        gnomad_afs = []
        for pop in ['gnomad_af_afr', 'gnomad_af_amr', 'gnomad_af_eas', 'gnomad_af_nfe', 'gnomad_af_sas']:
            if pop in variant and pd.notna(variant[pop]) and variant[pop] > 0:
                gnomad_afs.append(variant[pop])
        
        if gnomad_afs:
            return max(gnomad_afs)
        
        # Fallback to general gnomad_af
        if 'gnomad_af' in variant and pd.notna(variant['gnomad_af']):
            return variant['gnomad_af']
        
        return 0.0
        
    except Exception as e:
        logger.error(f"Error calculating max gnomAD AF: {e}")
        return 0.0

def format_cgen_frequency_with_counts(af, ac, an):
    """Format CGEN frequency with allele counts"""
    try:
        if pd.isna(af) or af == 0:
            af_str = "N/A"
        else:
            af_str = format_frequency(af)
        
        if pd.notna(ac) and pd.notna(an) and ac > 0 and an > 0:
            return f"{af_str} ({int(ac)} / {int(an)})"
        else:
            return af_str
            
    except Exception as e:
        logger.error(f"Error formatting CGEN frequency: {e}")
        return "N/A"

def validate_variant_data(variant_dict):
    """Validate variant data dictionary"""
    required_fields = ['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE']
    
    for field in required_fields:
        if field not in variant_dict or variant_dict[field] is None:
            return False, f"Missing required field: {field}"
    
    # Validate chromosome
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    chrom = str(variant_dict['CHROM']).replace('chr', '')
    if chrom not in valid_chroms:
        return False, f"Invalid chromosome: {chrom}"
    
    # Validate position
    try:
        pos = int(variant_dict['POS'])
        if pos <= 0:
            return False, "Position must be positive"
    except (ValueError, TypeError):
        return False, "Position must be a valid integer"
    
    # Validate alleles
    if not variant_dict['REF'] or not variant_dict['ALT']:
        return False, "REF and ALT alleles cannot be empty"
    
    return True, "Valid"

def create_variant_key(chrom, pos, ref, alt):
    """Create a unique variant key"""
    # Normalize chromosome name
    chrom = str(chrom).replace('chr', '')
    return f"{chrom}:{pos}:{ref}:{alt}"

def parse_variant_key(variant_key):
    """Parse variant key into components"""
    try:
        parts = variant_key.split(':')
        if len(parts) != 4:
            return None
        
        return {
            'CHROM': parts[0],
            'POS': int(parts[1]),
            'REF': parts[2],
            'ALT': parts[3]
        }
    except Exception:
        return None

def sort_variants(df, sort_by='position', ascending=True):
    """Sort variants by specified criteria"""
    if isinstance(df, pd.DataFrame):
        df = pl.from_pandas(df)
    
    if sort_by == 'position':
        # Convert chromosome to numeric for proper sorting
        df = df.with_columns([
            pl.when(pl.col('CHROM') == 'X').then(23)
            .when(pl.col('CHROM') == 'Y').then(24)
            .when(pl.col('CHROM') == 'MT').then(25)
            .otherwise(pl.col('CHROM').cast(pl.Int32, strict=False))
            .alias('chrom_numeric')
        ])
        df = df.sort(['chrom_numeric', 'POS'], descending=not ascending)
        df = df.drop('chrom_numeric')
    elif sort_by == 'gene':
        df = df.sort('gene', descending=not ascending)
    elif sort_by == 'consequence':
        df = df.sort('consequence', descending=not ascending)
    elif sort_by == 'clinvar':
        df = df.sort('clinvar_sig', descending=not ascending)
    elif sort_by == 'vaf':
        df = df.sort('VAF', descending=not ascending)
    elif sort_by == 'gnomad_af':
        # Sort by max_gnomad_af if available, otherwise gnomad_af
        if 'max_gnomad_af' in df.columns:
            df = df.sort('max_gnomad_af', descending=not ascending)
        else:
            df = df.sort('gnomad_af', descending=not ascending)
    
    return df

def get_severity_score(consequence):
    """Get numerical severity score for consequence"""
    severity_scores = {
        'transcript_ablation': 10,
        'splice_acceptor_variant': 9,
        'splice_donor_variant': 9,
        'stop_gained': 9,
        'stopgain': 9,
        'frameshift_variant': 9,
        'frameshift_deletion': 9,
        'frameshift_insertion': 9,
        'stop_lost': 8,
        'start_lost': 8,
        'transcript_amplification': 7,
        'inframe_insertion': 6,
        'inframe_deletion': 6,
        'missense_variant': 5,
        'nonsynonymous_SNV': 5,
        'protein_altering_variant': 5,
        'splice_region_variant': 4,
        'incomplete_terminal_codon_variant': 3,
        'start_retained_variant': 3,
        'stop_retained_variant': 3,
        'synonymous_variant': 2,
        'synonymous_SNV': 2,
        'coding_sequence_variant': 1,
        'variant': 0
    }
    
    return severity_scores.get(consequence, 0)

def get_clinvar_priority(significance):
    """Get numerical priority for ClinVar significance"""
    priority_scores = {
        'Pathogenic': 10,
        'Likely pathogenic': 9,
        'Drug response': 8,
        'Risk factor': 7,
        'Association': 6,
        'Protective': 5,
        'Conflicting': 4,
        'VUS': 3,
        'Uncertain significance': 3,
        'Likely benign': 2,
        'Benign': 1,
        'Other': 0,
        'Not provided': 0
    }
    
    return priority_scores.get(significance, 0)

def summarize_variants_by_gene(df):
    """Create gene-level summary of variants"""
    if isinstance(df, pd.DataFrame):
        df = pl.from_pandas(df)
    
    if len(df) == 0:
        return pl.DataFrame()
    
    gene_summary = (
        df.group_by('gene')
        .agg([
            pl.count('variant_key').alias('variant_count'),
            pl.col('SAMPLE').n_unique().alias('samples_affected'),
            pl.col('consequence').mode().first().alias('most_common_consequence'),
            pl.col('clinvar_sig').filter(pl.col('clinvar_sig').is_not_null()).mode().first().alias('most_common_clinvar'),
            pl.col('max_gnomad_af').max().alias('max_gnomad_af') if 'max_gnomad_af' in df.columns else pl.col('gnomad_af').max().alias('max_gnomad_af'),
            pl.col('cadd_score').max().alias('max_cadd_score')
        ])
        .sort('variant_count', descending=True)
    )
    
    return gene_summary

def create_summary_stats(df):
    """Create summary statistics for variants"""
    if isinstance(df, pd.DataFrame):
        df = pl.from_pandas(df)
    
    if len(df) == 0:
        return {}
    
    stats = {
        'total_variants': len(df),
        'unique_genes': df['gene'].n_unique(),
        'samples_with_variants': df['SAMPLE'].n_unique(),
        'consequence_distribution': df['consequence'].value_counts().to_pandas().to_dict(),
        'clinvar_distribution': df.filter(pl.col('clinvar_sig').is_not_null())['clinvar_sig'].value_counts().to_pandas().to_dict() if 'clinvar_sig' in df.columns else {},
        'chromosome_distribution': df['CHROM'].value_counts().to_pandas().to_dict(),
        'review_status_distribution': df['review_status'].value_counts().to_pandas().to_dict() if 'review_status' in df.columns else {}
    }
    
    return stats