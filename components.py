"""
UI Components for Variant Visualizer - MODIFIED VERSION
Modifications:
1. Hauteur fixe des lignes du tableau (48px)
2. Affichage compact des gènes dans le tableau (premier + nombre)
3. Affichage complet des gènes dans l'accordéon Clinical Information
"""

import dash_bootstrap_components as dbc
from dash import html, dcc
import pandas as pd
import polars as pl
from database import get_variant_comments, db
from utils import *
from config import *

# Import gene panel functions with fallback
try:
    from gene_panels import get_available_panels, get_panel_info
except ImportError:
    def get_available_panels():
        return []
    def get_panel_info(panel_id):
        return None

def create_gnomad_link(chrom, pos, ref, alt):
    """Create a link to gnomAD for a variant"""
    # gnomAD uses 1-based coordinates and specific format
    gnomad_url = f"https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}?dataset=gnomad_r2_1"
    
    return html.A(
        [html.I(className="fas fa-external-link-alt me-1"), "View in gnomAD"],
        href=gnomad_url,
        target="_blank",
        className="btn btn-outline-primary btn-sm",
        style={"fontSize": "12px"}
    )

def create_all_genes_display(gene_id_or_name):
    """Create display for all genes in Clinical Information section"""
    
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
        # Convert ID to name first
        gene_name = get_gene_name_from_id(gene)
        
        if gene_name == 'UNKNOWN' or not gene_name:
            # If we can't map it, show the original value
            gene_element = html.Span(
                gene, 
                style={
                    "color": "#6c757d", 
                    "fontStyle": "italic",
                    "fontSize": "14px",
                    "marginRight": "8px"
                }
            )
        else:
            # Create OMIM search URL using the gene name
            omim_url = f"https://www.omim.org/search?index=entry&start=1&limit=10&search={gene_name}"
            
            gene_element = html.A(
                gene_name,  # Always display the gene name, not the ID
                href=omim_url,
                target="_blank",
                style={
                    "fontWeight": "bold", 
                    "color": "#0097A7", 
                    "textDecoration": "none",
                    "fontSize": "14px",
                    "marginRight": "8px"
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
                    "fontSize": "14px",
                    "margin": "0 4px"
                })
            )
    
    # Return container with all genes
    return html.Div(
        gene_links,
        style={
            "display": "flex",
            "flexWrap": "wrap",
            "alignItems": "center",
            "gap": "2px",
            "lineHeight": "1.5"
        }
    )


def get_moi_badge(moi):
    """Return styled badge for Mode of Inheritance (MoI) - LARGER FONTS"""
    if not moi or pd.isna(moi) or moi in ['UNK', 'Unknown', '']:
        return dbc.Badge("UNK", color="secondary", 
                        style={"fontSize": "0.8em", "padding": "0.4em 0.6em"},
                        title="Mode of Inheritance Unknown")
    
    from config import MOI_COLORS
    color = MOI_COLORS.get(str(moi).upper(), 'secondary')
    
    # Tooltips explicatifs
    moi_tooltips = {
        'AD': 'Autosomal Dominant',
        'AR': 'Autosomal Recessive', 
        'XL': 'X-linked',
        'XLD': 'X-linked Dominant',
        'XLR': 'X-linked Recessive',
        'MT': 'Mitochondrial',
        'YL': 'Y-linked',
        'DD': 'Digenic',
        'OLI': 'Oligogenic',
        'SMU': 'Somatic Mutation',
        'UNK': 'Unknown'
    }
    
    tooltip = moi_tooltips.get(str(moi).upper(), f'Mode of Inheritance: {moi}')
    
    return dbc.Badge(
        str(moi).upper(), 
        color=color, 
        style={"fontSize": "0.8em", "padding": "0.4em 0.6em"},
        title=tooltip
    )

def create_aa_change_display(aa_change_str, variant_key, sample_id):
    """Create interactive display for AA changes with multiple transcripts and gene name mapping"""
    if not aa_change_str or aa_change_str in ['N/A', 'p.?', '.']:
        return html.Span("N/A", style={"fontSize": "14px", "color": "#6c757d"})
    
    # Split by comma to get multiple transcripts
    aa_changes = [change.strip() for change in str(aa_change_str).split(',') if change.strip()]

    # Extract gene ID from first part of each transcript and map to gene name
    aa_changes_mapped = []
    for change in aa_changes:
        parts = change.split(':', 1)  # split only at first colon
        gene_id = parts[0]
        rest = parts[1] if len(parts) > 1 else ''
        gene_name = get_gene_name_from_id(gene_id)
        aa_changes_mapped.append(f"{gene_name}:{rest}" if rest else gene_name)

    if len(aa_changes_mapped) <= 1:
        # Single transcript - display normally
        return html.Span(aa_changes_mapped[0], style={"fontSize": "14px"})
    
    # Multiple transcripts - create interactive display
    first_change = aa_changes_mapped[0]
    other_changes = aa_changes_mapped[1:]
    
    return html.Div([
        html.Span([
            html.Span(first_change, style={"fontSize": "14px", "fontWeight": "500"}),
            html.Span(
                f" (+{len(other_changes)} more)",
                style={
                    "fontSize": "12px", 
                    "color": "#0097A7", 
                    "fontWeight": "bold",
                    "marginLeft": "8px",
                    "cursor": "pointer",
                    "textDecoration": "underline"
                },
                title="Click to view other transcripts"
            )
        ], 
        id={"type": "aa-change-toggle", "variant": variant_key, "sample": sample_id},
        style={"cursor": "pointer"}),

        dbc.Collapse([
            dbc.Card([
                dbc.CardBody([
                    html.H6("Other Transcripts:", className="text-primary mb-2", style={"fontSize": "13px"}),
                    html.Div([
                        html.Div([
                            html.I(className="fas fa-dna me-2", style={"color": "#0097A7", "fontSize": "11px"}),
                            html.Span(change, style={"fontSize": "13px"})
                        ], className="mb-1")
                        for change in other_changes
                    ])
                ], style={"padding": "12px"})
            ], style={
                "marginTop": "8px", 
                "border": "1px solid #e9ecef",
                "borderRadius": "6px",
                "backgroundColor": "#f8f9fa"
            })
        ],
        id={"type": "aa-change-collapse", "variant": variant_key, "sample": sample_id},
        is_open=False)
    ])

def create_beautiful_variant_display(df):
    """Create optimized variant table display with FIXED ROW HEIGHT"""
    # Check DataFrame length properly for both Polars and Pandas
    if isinstance(df, pl.DataFrame):
        is_empty = len(df) == 0
        total_count = len(df)
    else:
        is_empty = df.empty
        total_count = len(df)
    
    if is_empty:
        return html.Div([
            dbc.Alert([
                html.Div([
                    html.I(className="fas fa-info-circle fa-3x text-info mb-3"),
                    html.H4("No Variants Found", className="text-info mb-3"),
                    html.P([
                        "No variants match your current selection or the database is empty. ",
                        html.Br(),
                        "Please check your filters or ensure your Parquet file contains data."
                    ], className="mb-0 text-muted")
                ], className="text-center")
            ], color="info", className="border-0", style={
                "background": "rgba(13, 202, 240, 0.1)",
                "borderRadius": "15px",
                "padding": "40px"
            })
        ])
    
    # Convert to pandas for display if it's Polars
    if isinstance(df, pl.DataFrame):
        df_pandas = df.to_pandas()
    else:
        df_pandas = df
    
    # OPTIMISATION MAJEURE: Limiter strictement à MAX_DISPLAY_VARIANTS
    display_df = df_pandas.head(MAX_DISPLAY_VARIANTS)
    showing_count = len(display_df)
    
    # Pre-calculate all data to avoid repeated calculations
    variant_data = []
    for idx, variant in display_df.iterrows():
        variant_key = variant.get('variant_key', f"{variant['CHROM']}:{variant['POS']}:{variant['REF']}:{variant['ALT']}")
        
        # Pre-format all data
        position = f"{variant['CHROM']}:{variant['POS']:}"
        variant_text = f"{variant['REF']}>{variant['ALT']}"
        vaf_display = format_percentage(variant.get('VAF', 0))
        
        # Use max_gnomad_af for display
        max_gnomad_af = variant.get('max_gnomad_af', variant.get('gnomad_af', 0))
        max_gnomad_af_display = format_frequency(max_gnomad_af)
        
        variant_data.append({
            'variant_key': variant_key,
            'position': position,
            'variant_text': variant_text,
            'vaf_display': vaf_display,
            'max_gnomad_af_display': max_gnomad_af_display,
            'max_gnomad_af': max_gnomad_af,
            'variant': variant
        })
    
    # Create table rows with FIXED HEIGHT - MODIFIED VERSION
    table_rows = []
    
    for data in variant_data:
        variant = data['variant']
        variant_key = data['variant_key']
        
        # Create data row with FIXED HEIGHT - 48px
        data_row = html.Tr([
            # Sample ID - FIXED HEIGHT
            html.Td(variant['SAMPLE'], style={
                "fontSize": "14px", "fontWeight": "600", "height": "48px", 
                "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Position - FIXED HEIGHT
            html.Td(data['position'], style={
                "fontFamily": "monospace", "fontSize": "14px", "fontWeight": "500",
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Gene - FIXED HEIGHT AND SINGLE LINE (COMPACT VERSION)
            html.Td(
                html.Div(create_gene_link(variant.get('gene', 'UNKNOWN')), style={
                    "whiteSpace": "nowrap", "overflow": "hidden", "textOverflow": "ellipsis",
                    "maxWidth": "120px"
                }), 
                style={
                    "fontSize": "14px", "height": "48px", "verticalAlign": "middle", 
                    "padding": "8px 15px"
                }
            ),
            # Variant - FIXED HEIGHT
            html.Td(data['variant_text'], style={
                "fontFamily": "monospace", "fontSize": "14px", "backgroundColor": "#f8f9fa", 
                "borderRadius": "4px", "textAlign": "center", "fontWeight": "500",
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Genotype - FIXED HEIGHT
            html.Td([get_genotype_badge_optimized(variant.get('GT', './.'))], style={
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # NOUVEAU : MoI - FIXED HEIGHT
            html.Td([get_moi_badge(variant.get('moi', 'UNK'))], style={
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # VAF - FIXED HEIGHT
            html.Td(data['vaf_display'], style={
                "fontFamily": "monospace", "fontSize": "13px", "textAlign": "right", 
                "fontWeight": "500", "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Consequence - FIXED HEIGHT
            html.Td([get_consequence_badge_optimized(variant.get('consequence', 'variant'))], style={
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # ClinVar Classification - FIXED HEIGHT
            html.Td([get_clinvar_badge_optimized(variant.get('clinvar_sig'))], style={
                "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Population Frequency (gnomAD max) - FIXED HEIGHT
            html.Td([
                html.Span(data['max_gnomad_af_display'], 
                    style={
                        "fontFamily": "monospace", 
                        "fontSize": "13px", 
                        "fontWeight": "bold",
                        "color": "#dc3545" if data['max_gnomad_af'] and data['max_gnomad_af'] < 0.001 
                            else "#ffc107" if data['max_gnomad_af'] and 0.001 <= data['max_gnomad_af'] < 0.01 
                            else "#28a745" if data['max_gnomad_af'] and data['max_gnomad_af'] >= 0.01 
                            else "#6c757d"
                    },
                    title=f"Max gnomAD Population Allele Frequency: {data['max_gnomad_af_display']}"
                )
            ], style={
                "textAlign": "right", "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
            # Comments - FIXED HEIGHT
            html.Td([
                dbc.Button([html.I(className="fas fa-comments me-1"), str(variant.get('comment_count', 0))], 
                        id={"type": "comment-btn", "variant": variant_key, "sample": variant['SAMPLE']},
                        color="outline-primary", size="sm", n_clicks=0, 
                        style={"fontSize": "12px"})
            ], style={
                "textAlign": "center", "height": "48px", "verticalAlign": "middle", "padding": "8px 15px"
            }),
        ], className="variant-row", id={"type": "variant-row", "variant": variant_key, "sample": variant['SAMPLE']}, n_clicks=0)
        
        # Accordion collapse row - LAZY LOADING
        accordion_row = html.Tr([
            html.Td([
                dbc.Collapse([
                    html.Div(id={"type": "variant-details-lazy", "variant": variant_key, "sample": variant['SAMPLE']})
                ], 
                id={"type": "variant-collapse", "variant": variant_key, "sample": variant['SAMPLE']},
                is_open=False)
            ], colSpan=10, style={"padding": "0", "border": "none"})
        ], style={"border": "none"})
        
        table_rows.extend([data_row, accordion_row])
    
    # Performance notice
    performance_notice = []
    if total_count > MAX_DISPLAY_VARIANTS:
        performance_notice = [
            dbc.Alert([
                html.I(className="fa fa-exclamation-triangle"),
                f" Performance mode: Showing first {showing_count} of {total_count:,} variants. ",
                
            ], color="warning", className="mb-3", style={"fontSize": "0.9em"})
        ]
    
    # Create table with minimal overhead - LARGER HEADERS
    return html.Div([
        *performance_notice,
        html.Table([
            html.Thead([
                html.Tr([
                    html.Th("Sample", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Position", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Gene", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Variant", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Genotype", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
					html.Th("MoI", className="sortable-header", title="Mode of Inheritance", style={"fontSize": "15px", "fontWeight": "700"}),  # NOUVEAU
                    html.Th("VAF", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Consequence", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("ClinVar", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("AF (gnomAD)", className="sortable-header", title="Max Population Allele Frequency from gnomAD", style={"fontSize": "15px", "fontWeight": "700"}),
                    html.Th("Comments", className="sortable-header", style={"fontSize": "15px", "fontWeight": "700"}),
                ])
            ]),
            html.Tbody(table_rows)
        ], className="variants-table")
    ], className="variants-table-container")

def create_variant_details_accordion(variant):
    """Create detailed information panel for accordion expansion with ALL GENES in Clinical Information"""
    
    variant_key = variant.get('variant_key', f"{variant['CHROM']}:{variant['POS']}:{variant['REF']}:{variant['ALT']}")
    sample_id = variant['SAMPLE']
    
    # Get comments for this variant and sample
    try:
        comments_df = get_variant_comments(variant_key, sample_id)
    except Exception as e:
        logger.error(f"Error retrieving comments: {e}")
        comments_df = pd.DataFrame()
    
    return html.Div([
        dbc.Row([
            # Column 1: Variant Information
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-dna me-2"),
                        "Variant Details"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        html.Div([
                            html.Strong("Position: "),
                            html.Span(f"{variant['CHROM']}:{variant['POS']:}", style={"fontSize": "14px"})
                        ], className="mb-2"),

                        html.Div([
                            html.Strong("Consequence: "),
                            get_consequence_badge(variant.get('consequence', 'variant'))
                        ], className="mb-2"),

                        html.Div([
                            html.Strong("Genotype: "),
                            get_genotype_badge(variant.get('GT', './.'))
                        ], className="mb-2"),

                        html.Div([
                            html.Strong("AA Change: "),
                            create_aa_change_display(
                                variant.get('aa_change', 'N/A'), 
                                variant_key, 
                                sample_id
                            )
                        ], className="mb-2"),

                    ])
                ], className="detail-section uniform-height")
            ], width=4),
            
            # Column 2: Technical Information
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-microscope me-2"),
                        "Technical Details"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        html.Div([
                            html.Strong("VAF: "),
                            html.Span(format_percentage(variant.get('VAF', 0)), style={"fontSize": "14px"})
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("Depth: "),
                            html.Span(f"{variant.get('DP', 0)}" if pd.notna(variant.get('DP')) else "N/A", style={"fontSize": "14px"})
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("Allelic Depth: "),
                            html.Span(variant.get('AD', 'N/A'), style={"fontSize": "14px"})
                        ], className="mb-2"),

                        html.Div([
                            html.Strong("Quality: "),
                            html.Span(f"{variant.get('QUAL', 0):.1f}" if pd.notna(variant.get('QUAL')) else "N/A", style={"fontSize": "14px"})
                        ], className="mb-2"),
                    ])
                ], className="detail-section uniform-height")
            ], width=4),
            
            # Column 3: Clinical Information - MODIFIED TO SHOW ALL GENES AND MoI
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-stethoscope me-2"),
                        "Clinical Information"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        html.Div([
                            html.Strong("ClinVar Significance: "),
                            get_clinvar_badge(variant.get('clinvar_sig'))
                        ], className="mb-2"),
                        
                        # Genes (existant)
                        html.Div([
                            html.Strong("Gene(s): ", style={"marginRight": "8px"}),
                            create_all_genes_display(variant.get('gene', 'UNKNOWN'))
                        ], className="mb-2", style={"display": "flex", "alignItems": "flex-start", "flexWrap": "wrap"}),

                        # NOUVEAU : Mode d'héritage
                        html.Div([
                            html.Strong("Mode of Inheritance: "),
                            get_moi_badge(variant.get('moi', 'UNK'))
                        ], className="mb-2"),

                    ])
                ], className="detail-section uniform-height")
            ], width=4)
        ]),  # CORRECTION: Fermeture du premier dbc.Row
        
        # Bioinformatics Scores Section
        dbc.Row([
            # Column 1: Pathogenicity Prediction
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-calculator me-2"),
                        "Pathogenicity Prediction"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        html.Div([
                            html.Strong("CADD Score: "),
                            html.Span(
                                format_score(variant.get('cadd_score', 0)) if variant.get('cadd_score') else "N/A",
                                className="fw-bold",
                                style={"color": "#28a745" if variant.get('cadd_score', 0) and variant.get('cadd_score', 0) < 15 
                                    else "#ffc107" if variant.get('cadd_score', 0) and 15 <= variant.get('cadd_score', 0) < 25 
                                    else "#dc3545" if variant.get('cadd_score', 0) and variant.get('cadd_score', 0) >= 25 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (>15 deleterious)", className="text-muted")
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("SIFT Score: "),
                            html.Span(
                                format_score(variant.get('sift_score', 0)) if variant.get('sift_score') else "N/A",
                                className="fw-bold",
                                style={"color": "#dc3545" if variant.get('sift_score', 0) and variant.get('sift_score', 0) < 0.05 
                                    else "#28a745" if variant.get('sift_score', 0) and variant.get('sift_score', 0) >= 0.05 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (<0.05 deleterious)", className="text-muted")
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("PolyPhen Score: "),
                            html.Span(
                                format_score(variant.get('polyphen_score', 0)) if variant.get('polyphen_score') else "N/A",
                                className="fw-bold",
                                style={"color": "#28a745" if variant.get('polyphen_score', 0) and variant.get('polyphen_score', 0) < 0.15 
                                    else "#ffc107" if variant.get('polyphen_score', 0) and 0.15 <= variant.get('polyphen_score', 0) < 0.85 
                                    else "#dc3545" if variant.get('polyphen_score', 0) and variant.get('polyphen_score', 0) >= 0.85 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (>0.85 damaging)", className="text-muted")
                        ], className="mb-2")
                    ])
                ], className="detail-section uniform-height")
            ], width=4),
            
            # Column 2: Splicing & Conservation Scores
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-dna me-2"),
                        "Splicing & Conservation"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        html.Div([
                            html.Strong("REVEL Score: "),
                            html.Span(
                                format_score(variant.get('revel_score', 0)) if variant.get('revel_score') else "N/A",
                                className="fw-bold",
                                style={"color": "#28a745" if variant.get('revel_score', 0) and variant.get('revel_score', 0) < 0.3 
                                    else "#ffc107" if variant.get('revel_score', 0) and 0.3 <= variant.get('revel_score', 0) < 0.7 
                                    else "#dc3545" if variant.get('revel_score', 0) and variant.get('revel_score', 0) >= 0.7 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (>0.7 pathogenic)", className="text-muted")
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("SpliceAI Score: "),
                            html.Span(
                                format_score(variant.get('splice_ai', 0)) if variant.get('splice_ai') else "N/A",
                                className="fw-bold",
                                style={"color": "#28a745" if variant.get('splice_ai', 0) and variant.get('splice_ai', 0) < 0.2 
                                    else "#ffc107" if variant.get('splice_ai', 0) and 0.2 <= variant.get('splice_ai', 0) < 0.5 
                                    else "#dc3545" if variant.get('splice_ai', 0) and variant.get('splice_ai', 0) >= 0.5 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (>0.5 splice altering)", className="text-muted")
                        ], className="mb-2"),
                        
                        html.Div([
                            html.Strong("pLI Score: "),
                            html.Span(
                                format_score(variant.get('pli_score', 0)) if variant.get('pli_score') else "N/A",
                                className="fw-bold",
                                style={"color": "#28a745" if variant.get('pli_score', 0) and variant.get('pli_score', 0) < 0.1 
                                    else "#ffc107" if variant.get('pli_score', 0) and 0.1 <= variant.get('pli_score', 0) < 0.9 
                                    else "#dc3545" if variant.get('pli_score', 0) and variant.get('pli_score', 0) >= 0.9 
                                    else "#6c757d", "fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Small(" (>0.9 intolerant)", className="text-muted")
                        ], className="mb-2")
                    ])
                ], className="detail-section uniform-height")
            ], width=4),
            
            # Column 3: Population & Additional Scores
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-globe me-2"),
                        "Population Data"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        # Display max gnomAD AF
                        html.Div([
                            html.Strong("Allele Frequency gnomAD (max): "),
                            html.Span(
                                format_frequency(variant.get('max_gnomad_af', 0)) if variant.get('max_gnomad_af') else "N/A",
                                className="fw-bold",
                                style={"color": "#dc3545" if variant.get('max_gnomad_af', 0) and variant.get('max_gnomad_af', 0) < 0.001 
                                    else "#ffc107" if variant.get('max_gnomad_af', 0) and 0.001 <= variant.get('max_gnomad_af', 0) < 0.01 
                                    else "#28a745" if variant.get('max_gnomad_af', 0) and variant.get('max_gnomad_af', 0) >= 0.01 
                                    else "#6c757d", "fontSize": "14px"}
                            )
                        ], className="mb-2"),
                        
                        # Allele Count gnomAD with homo/hemi counts on separate lines
                        html.Div([
                            html.Strong("Allele Count gnomAD: "),
                            html.Span(
                                f"{variant.get('ac_gnomad', 0)}" if pd.notna(variant.get('ac_gnomad')) else "N/A",
                                style={"fontSize": "14px"}
                            ),
                            html.Br(),
                            html.Span([
                                html.Span("             • ", style={"fontSize": "14px"}),
                                html.Strong("NB of homo : "),
                                f"{variant.get('nhomalt_gnomad', 0):.0f}" if pd.notna(variant.get('nhomalt_gnomad')) else "N/A"
                            ], style={"fontSize": "14px"}),
                            html.Br(),
                            html.Span([
                                html.Span("             • ", style={"fontSize": "14px"}),
                                html.Strong("NB of hemi : "),
                                f"{variant.get('nhemalt_gnomad', 0):.0f}" if pd.notna(variant.get('nhemalt_gnomad')) else "N/A"
                            ], style={"fontSize": "14px"})
                        ], className="mb-2"),
                        
                        # AF_CGEN with calculation and color coding (only frequency colored)
                        html.Div([
                            html.Strong("Allele Frequency CGEN: "),
                            html.Span([
                                html.Span(
                                    f"{format_frequency(variant.get('af_cgen', 0))}" if variant.get('af_cgen') else "N/A",
                                    className="fw-bold",
                                    style={"color": "#dc3545" if variant.get('af_cgen', 0) and variant.get('af_cgen', 0) < 0.001 
                                        else "#ffc107" if variant.get('af_cgen', 0) and 0.001 <= variant.get('af_cgen', 0) < 0.01 
                                        else "#28a745" if variant.get('af_cgen', 0) and variant.get('af_cgen', 0) >= 0.01 
                                        else "#6c757d", "fontSize": "14px"}
                                ),
                                html.Span(
                                    f" ({variant.get('ac_cgen', 0)}/{variant.get('an_cgen', 0)})" if variant.get('ac_cgen') and variant.get('an_cgen') else "",
                                    style={"fontSize": "14px", "color": "#6c757d"}
                                )
                            ])
                        ], className="mb-3"),
                        
                        # gnomAD link
                        html.Div([
                            create_gnomad_link(variant['CHROM'], variant['POS'], variant['REF'], variant['ALT'])
                        ])
                    ])
                ], className="detail-section uniform-height")
            ], width=4)
        ]),  # CORRECTION: Fermeture du deuxième dbc.Row
            
        # Comments Section
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H6([
                        html.I(className="fas fa-comments me-2"),
                        f"Comments ({len(comments_df)})"
                    ], className="text-primary mb-3"),
                    
                    html.Div([
                        create_comments_display_accordion(comments_df, variant_key, sample_id)
                    ], className="comments-section")
                ], className="detail-section")
            ], width=12)
        ])  # CORRECTION: Fermeture du troisième dbc.Row
    ], style={"padding": "20px", "backgroundColor": "rgba(240, 248, 255, 0.8)", "borderRadius": "0 0 12px 12px"})

def create_database_status_display():
	"""Create database status display with optimized queries"""
	stats = db.get_database_stats()
	
	if stats['total_variants'] == 0:
		return html.Div([
			dbc.Row([
				dbc.Col([html.I(className="fas fa-info-circle me-2 text-info"), html.Strong("0"), html.Span(" variants", className="ms-1")], width="auto"),
				dbc.Col([html.I(className="fas fa-users me-2 text-secondary"), html.Strong("0"), html.Span(" samples", className="ms-1")], width="auto"),
				dbc.Col([dbc.Badge("Empty Database", color="light", className="ms-2")], width="auto"),
				dbc.Col([dbc.Badge("Parquet DB", color="success", className="ms-2")], width="auto")
			], align="center", className="g-3")
		])
	
	return html.Div([
		dbc.Row([
			dbc.Col([html.I(className="fas fa-table me-2 text-primary"), html.Strong(f"{stats['total_variants']:,}"), html.Span(" variants", className="ms-1")], width="auto"),
			dbc.Col([html.I(className="fas fa-users me-2 text-info"), html.Strong(f"{stats['total_samples']}"), html.Span(" samples", className="ms-1")], width="auto"),
			dbc.Col([html.I(className="fas fa-dna me-2 text-success"), html.Strong(f"{stats['clinvar_annotated']}"), html.Span(" ClinVar", className="ms-1")], width="auto"),
			dbc.Col([html.I(className="fas fa-check-circle me-2 text-success"), html.Strong(f"{stats['reviewed_variants']}"), html.Span(" reviewed", className="ms-1")], width="auto"),
			dbc.Col([html.I(className="fas fa-clock me-2 text-warning"), html.Strong(f"{stats['pending_variants']}"), html.Span(" pending", className="ms-1")], width="auto"),
			dbc.Col([dbc.Badge("Parquet DB", color="success", className="ms-2")], width="auto")
		], align="center", className="g-3")
	])

def create_gene_panel_selector():
	"""Create the gene panel selector component - SIMPLIFIED (no green gene filter) - HARMONIZED STYLE"""
	available_panels = get_available_panels()
	
	return dbc.Card([
		dbc.CardBody([
			html.Div([
				html.Div([
					html.H6([html.I(className="fas fa-list-alt me-2"), "Gene Panel Selection"], 
						   className="mb-2 text-primary", style={"fontSize": "1.1rem"}),
					html.P("Select gene panels to filter variants by specific gene sets:", 
						  className="text-muted mb-0", style={"fontSize": "14px"})
				]),
				html.Div([
					dcc.Dropdown(
						id="gene-panel-selector",
						options=available_panels,
						placeholder="Select gene panels...",
						multi=True,
						searchable=True,
						clearable=True,
						style={"minWidth": "400px", "fontSize": "15px", "zIndex": "9999"}
					)
				], style={"flex": "1", "position": "relative", "zIndex": "9999"}),
				dbc.ButtonGroup([
					dbc.Button([
						html.I(className="fas fa-info-circle me-1"),
						"Panel Info"
					], id="panel-info-btn", color="outline-info", size="sm",
					style={"fontSize": "13px"}),
					dbc.Button([
						html.I(className="fas fa-sync me-1"),
						"Update Panels"
					], id="update-panels-btn", color="outline-secondary", size="sm",
					style={"fontSize": "13px"}),
					dbc.Button("Clear", id="clear-gene-panels", color="outline-secondary", size="sm",
							  style={"fontSize": "13px"})
				])
			], style={"display": "flex", "alignItems": "center", "gap": "15px"})
		], style={"padding": "1.5rem"})
	], className="glass-card mb-3", style={"position": "relative", "zIndex": "1000"})

def create_panel_info_modal():
	"""Create modal for displaying panel information"""
	return dbc.Modal([
		dbc.ModalHeader([
			dbc.ModalTitle(id="panel-info-modal-title", children="Panel Information")
		]),
		dbc.ModalBody([
			html.Div(id="panel-info-content")
		]),
		dbc.ModalFooter([
			dbc.Button("Close", id="close-panel-info-modal", className="ms-auto", color="secondary")
		])
	], id="panel-info-modal", is_open=False, size="lg")

def create_panel_info_content(panel_ids, selected_panels):
	"""Create content for panel info modal - SIMPLIFIED (no green gene specific info)"""
	if not panel_ids:
		return html.Div([
			dbc.Alert("No panels selected.", color="info")
		])
	
	panel_info_cards = []
	total_genes = set()
	
	for panel_id in panel_ids:
		panel_info = get_panel_info(panel_id)
		if panel_info:
			# Get genes for this panel
			try:
				from gene_panels import get_genes_for_panels
				all_panel_genes = get_genes_for_panels([panel_id])
				total_genes.update(all_panel_genes)
			except ImportError:
				all_panel_genes = []
			
			source_icon = {
				'panelapp_uk': '🇬🇧',
				'panelapp_au': '🇦🇺', 
				'internal': '🇨🇭'
			}.get(panel_info.get('source', ''), '📋')
			
			# Show gene confidence info for external panels
			gene_confidence_info = html.Div()
			if panel_info.get('source') in ['panelapp_uk', 'panelapp_au']:
				green_count = panel_info.get('green_gene_count', 0)
				amber_count = panel_info.get('amber_gene_count', 0)
				red_count = panel_info.get('red_gene_count', 0)
				unknown_count = panel_info.get('gene_count', 0) - green_count - amber_count - red_count
				
				badges = []
				if green_count > 0:
					badges.append(dbc.Badge(f"{green_count} Green", color="success", className="me-1"))
				if amber_count > 0:
					badges.append(dbc.Badge(f"{amber_count} Amber", color="warning", className="me-1"))
				if red_count > 0:
					badges.append(dbc.Badge(f"{red_count} Red", color="danger", className="me-1"))
				if unknown_count > 0:
					badges.append(dbc.Badge(f"{unknown_count} Unknown", color="secondary", className="me-1"))
				
				gene_confidence_info = html.Div([
					html.Div([
						html.Strong("Gene Confidence Distribution: "),
						html.Div(badges, style={"display": "inline-block"})
					], className="mb-2")
				])
			
			card = dbc.Card([
				dbc.CardHeader([
					html.H6([
						source_icon, " ", panel_info.get('panel_name', 'Unknown Panel')
					], className="mb-0 text-primary")
				]),
				dbc.CardBody([
					html.Div([
						html.Strong("Source: "), 
						dbc.Badge(panel_info.get('source', 'Unknown'), color="info", className="me-2")
					], className="mb-2"),
					html.Div([
						html.Strong("Version: "), 
						panel_info.get('panel_version', 'Unknown')
					], className="mb-2"),
					gene_confidence_info,
					html.Div([
						html.A("View Online", 
							  href=panel_info.get('panel_url', '#'), 
							  target="_blank",
							  className="btn btn-sm btn-outline-primary") if panel_info.get('panel_url') else html.Div()
					])
				])
			], className="mb-3")
			
			panel_info_cards.append(card)
	
	# Summary with enhanced confidence information - COUNT UNIQUE GENES CORRECTLY
	all_genes_with_confidence = {}  # gene_symbol -> confidence_level
	
	for panel_id in panel_ids:
		panel_info = get_panel_info(panel_id)
		if panel_info:
			# Get all genes for this panel to collect unique genes with their confidence
			try:
				from gene_panels import panel_manager
				panel_genes = panel_manager.panels_df.filter(pl.col('panel_id') == panel_id)
				
				# For each gene, store its confidence (latest panel wins if conflict)
				for _, row in panel_genes.to_pandas().iterrows():
					gene_symbol = row['gene_symbol']
					confidence = row['gene_confidence']
					
					# If gene already exists, prioritize: GREEN > AMBER > RED > UNKNOWN
					if gene_symbol in all_genes_with_confidence:
						current_conf = all_genes_with_confidence[gene_symbol]
						priority = {'GREEN': 4, 'AMBER': 3, 'RED': 2, 'UNKNOWN': 1}
						if priority.get(confidence, 0) > priority.get(current_conf, 0):
							all_genes_with_confidence[gene_symbol] = confidence
					else:
						all_genes_with_confidence[gene_symbol] = confidence
						
			except Exception as e:
				logger.error(f"Error getting unique genes for panel {panel_id}: {e}")
	
	# Count unique genes by confidence
	unique_green_genes = set()
	unique_amber_genes = set()
	unique_red_genes = set()
	unique_unknown_genes = set()
	
	for gene_symbol, confidence in all_genes_with_confidence.items():
		if confidence == 'GREEN':
			unique_green_genes.add(gene_symbol)
		elif confidence == 'AMBER':
			unique_amber_genes.add(gene_symbol)
		elif confidence == 'RED':
			unique_red_genes.add(gene_symbol)
		elif confidence == 'UNKNOWN':
			unique_unknown_genes.add(gene_symbol)
	
	total_unique_green = len(unique_green_genes)
	total_unique_amber = len(unique_amber_genes)
	total_unique_red = len(unique_red_genes)
	total_unique_unknown = len(unique_unknown_genes)
	
	# Verify the math: total unique genes from confidence should match total_genes
	calculated_total = total_unique_green + total_unique_amber + total_unique_red + total_unique_unknown
	logger.info(f"Gene confidence math check: {total_unique_green} + {total_unique_amber} + {total_unique_red} + {total_unique_unknown} = {calculated_total}, expected: {len(total_genes)}")
	
	# Enhanced summary with confidence distribution - ONLY SHOW IF WE HAVE EXTERNAL PANELS
	confidence_summary = html.Div()
	has_external_panels = any(panel_info.get('source') in ['panelapp_uk', 'panelapp_au'] for panel_id in panel_ids for panel_info in [get_panel_info(panel_id)] if panel_info)
	
	if has_external_panels and (total_unique_green > 0 or total_unique_amber > 0 or total_unique_red > 0):
		confidence_badges = []
		if total_unique_green > 0:
			confidence_badges.append(dbc.Badge(f"{total_unique_green} Green genes", color="success", className="me-1"))
		if total_unique_amber > 0:
			confidence_badges.append(dbc.Badge(f"{total_unique_amber} Amber genes", color="warning", className="me-1"))
		if total_unique_red > 0:
			confidence_badges.append(dbc.Badge(f"{total_unique_red} Red genes", color="danger", className="me-1"))
		if total_unique_unknown > 0:
			confidence_badges.append(dbc.Badge(f"{total_unique_unknown} Unknown genes", color="secondary", className="me-1"))
		
		if confidence_badges:  # Only show if we have badges to show
			confidence_summary = html.Div([
				html.Div([
					html.Strong("Overall Gene Confidence Summary: "),
					html.Div(confidence_badges, style={"display": "inline-block"})
				], className="mb-2"),
				html.Small([
					"",
					#html.Span(f"Total: {calculated_total} genes", className="fw-bold")
				], className="text-muted d-block", style={"fontSize": "12px"})
			])
	elif has_external_panels:
		# Show debug info if we have external panels but no confidence data
		confidence_summary = html.Div([
			html.Div([
				html.Strong("Gene Confidence Summary: "),
				dbc.Badge(f"All genes marked as Unknown - debugging in progress", color="warning")
			], className="mb-2"),
			html.Small("Check console logs for confidence parsing debug information", className="text-muted d-block", style={"fontSize": "12px"})
		])
	
	summary_card = dbc.Card([
		dbc.CardHeader([
			html.H6("Summary", className="mb-0 text-success")
		]),
		dbc.CardBody([
			html.Div([
				html.Strong("Total Panels Selected: "), f"{len(panel_ids)}"
			], className="mb-2"),
			html.Div([
				html.Strong("Total Unique Genes: "), f"{len(total_genes)}"
			], className="mb-2"),
			confidence_summary,
			html.Details([
				html.Summary("Show All Genes", style={"cursor": "pointer", "fontWeight": "bold"}),
				html.Div([
					html.Div([
						dbc.Badge(gene, color="info", className="me-1 mb-1") 
						for gene in sorted(total_genes)
					])
				], style={"maxHeight": "300px", "overflowY": "auto", "marginTop": "10px"})
			]) if total_genes else html.Div()
		])
	], className="mb-3", color="success", outline=True)
	
	return html.Div([summary_card] + panel_info_cards)

def create_update_status_toast():
	"""Create toast notification for panel updates"""
	return dbc.Toast(
		id="update-status-toast",
		header="Panel Update Status",
		is_open=False,
		dismissable=True,
		duration=4000,
		icon="info",
		style={"position": "fixed", "top": 66, "right": 10, "width": 350, "zIndex": 9999}
	)

def create_comments_display_accordion(comments_df, variant_key, sample_id):
	"""Create comments display section for accordion"""
	if comments_df.empty:
		return html.Div([
			html.P("No comments yet.", className="text-muted text-center"),
			dbc.Button([
				html.I(className="fas fa-plus me-2"),
				"Add First Comment"
			], 
			id={"type": "add-comment-btn", "variant": variant_key, "sample": sample_id},
			color="primary",
			size="sm",
			outline=True,
			className="d-block mx-auto")
		])
	
	comment_items = []
	for _, comment in comments_df.iterrows():
		comment_item = html.Div([
			html.Div([
				html.Strong([
					html.I(className="fas fa-user-circle me-2"),
					comment['user_name']
				], className="text-primary"),
				html.Small(comment['timestamp'], className="text-muted float-end")
			], className="mb-2"),
			html.P(comment['comment_text'], className="mb-0")
		], className="comment-item")
		comment_items.append(comment_item)
	
	comment_items.append(
		dbc.Button([
			html.I(className="fas fa-plus me-2"),
			"Add Comment"
		], 
		id={"type": "add-comment-btn", "variant": variant_key, "sample": sample_id},
		color="primary",
		size="sm",
		outline=True,
		className="mt-2")
	)
	
	return html.Div(comment_items)

def create_sidebar():
	"""Create the advanced filters sidebar with Search, VAF Range, and Apply button"""
	return html.Div([
		# Sidebar header
		html.Div([
			html.H4([
				html.I(className="fas fa-sliders-h me-2"),
				"Advanced Filters"
			], className="text-primary mb-0", style={"fontSize": "1.3rem"}),
			dbc.Button([
				html.I(className="fas fa-times")
			], 
			id="close-sidebar-btn",
			color="light",
			size="sm",
			className="ms-auto")
		], className="d-flex align-items-center justify-content-between p-3 border-bottom"),
		
		# Sidebar content
		html.Div([
			# Search
			html.Div([
				html.Label([
					html.I(className="fas fa-search me-2"),
					"Search Genes, Samples..."
				], className="fw-bold mb-2 text-primary", style={"fontSize": "15px"}),
				dbc.Input(
					id="search-input",
					placeholder="e.g., BRCA1, Sample001, chr1...",
					debounce=True,
					style={"borderRadius": "8px", "fontSize": "15px", "padding": "0.6rem"}
				),
				html.Small(
					"Search by gene name, sample ID, chromosome, or consequence", 
					className="text-muted mt-1 d-block", style={"fontSize": "13px"}
				)
			], className="mb-4"),
			
			# VAF Range Filter
			html.Div([
				html.Label([
					html.I(className="fas fa-chart-line me-2"),
					"Variant Allele Frequency (VAF)"
				], className="fw-bold mb-2 text-primary", style={"fontSize": "15px"}),
				html.Label("Filter by VAF Range:", className="mb-2", 
						  style={"fontSize": "14px"}),
				dcc.RangeSlider(
					id="vaf-range-slider",
					min=0,
					max=1,
					step=0.01,
					value=[0, 1],
					marks={0: '0%', 0.25: '25%', 0.5: '50%', 0.75: '75%', 1: '100%'},
					tooltip={"placement": "bottom", "always_visible": True}
				),
				html.Small(
					"Only show variants within the selected VAF range", 
					className="text-muted mt-2 d-block", style={"fontSize": "13px"}
				)
			], className="mb-4"),
			
			# Apply and Clear buttons
			html.Div([
				dbc.Button([
					html.I(className="fas fa-search me-2"),
					"Apply Filters & Search"
				], 
				id="apply-filters-btn", 
				color="primary", 
				size="md", 
				className="w-100 mb-3",
				style={"borderRadius": "8px", "fontWeight": "bold", "fontSize": "14px", "padding": "0.7rem"}
				),
				dbc.Button([
					html.I(className="fas fa-trash me-2"),
					"Clear All Filters"
				], 
				id="clear-filters", 
				color="outline-danger", 
				size="sm", 
				className="w-100",
				style={"borderRadius": "8px", "fontSize": "13px"}
				)
			])
		], className="p-3")
	], id="filter-sidebar", className="filter-sidebar")

def create_header():
	"""Create the main header component"""
	return dbc.Card([
		dbc.CardBody([
			html.Div([
				html.Div([
					html.H1([
						"🧬 ", 
						html.Span("Variant Visualizer", style={"color": "#0097A7"}),
					], className="mb-0", style={"fontSize": "2.2rem"}),
				], style={"flex": "1"}),
				html.Div([
					# Space for future buttons if needed
				])
			], style={"display": "flex", "alignItems": "center", "justifyContent": "space-between"})
		], style={"padding": "1.5rem"})
	], className="glass-card mb-3")

def create_sample_selector():
	"""Create the sample selector component"""
	return dbc.Card([
		dbc.CardBody([
			html.Div([
				html.Div([
					html.H6([html.I(className="fas fa-users me-2"), "Sample Selection"], 
						   className="mb-2 text-primary", style={"fontSize": "1.1rem"}),
					html.P("Select samples to display variants:", 
						  className="text-muted mb-0", style={"fontSize": "14px"})
				]),
				html.Div([
					dcc.Dropdown(
						id="sample-selector",
						placeholder="Select samples...",
						multi=True,
						searchable=True,
						clearable=True,
						style={"minWidth": "300px", "fontSize": "15px", "zIndex": "9999"}
					)
				], style={"flex": "1", "position": "relative", "zIndex": "9999"}),
				dbc.ButtonGroup([
					dbc.Button("Clear", id="clear-samples", color="outline-secondary", size="sm",
							  style={"fontSize": "13px"})
				])
			], style={"display": "flex", "alignItems": "center", "gap": "15px"})
		], style={"padding": "1.5rem"})
	], className="glass-card mb-3", style={"position": "relative", "zIndex": "1000"})

def create_main_filters_panel():
	"""Create the main filters panel with quick filters and reset button"""
	return dbc.Card([
		dbc.CardBody([
			dbc.Row([
				dbc.Col([
					html.Div(id="variant-count")
				], width=3),
				dbc.Col([
					html.Div([
						html.Label("", className="fw-bold me-3 mb-0",
								  style={"fontSize": "15px"}),
						html.Div([
							dbc.Button([f["icon"], " ", f["name"]], 
									id={"type": "preset-filter", "id": f["id"]},
									color="outline-info", size="sm", className="quick-filter-btn me-2 mb-1", 
									n_clicks=0, style={"fontSize": "14px"})
							for f in PRESET_FILTERS
						], style={"display": "flex", "flexWrap": "wrap"})
					])
				], width=6),
				dbc.Col([
					dbc.ButtonGroup([
						dbc.Button("More Filters", id="more-filters-btn", color="outline-secondary", size="sm",
								  style={"fontSize": "13px"}),
						dbc.Button([
							html.I(className="fas fa-undo me-1"),
							"Reset"
						], id="reset-all-btn", color="outline-danger", size="sm", title="Reset all filters",
						style={"fontSize": "13px"})
					])
				], width=3, className="text-end")
			], align="center")
		], style={"padding": "1.25rem"})
	], className="glass-card mb-3")

def create_comment_modal():
	"""Create the comment modal component"""
	return dbc.Modal([
		dbc.ModalHeader([
			dbc.ModalTitle(id="comment-modal-title")
		]),
		dbc.ModalBody([
			html.Div(id="existing-comments"),
			html.Hr(),
			html.H6("💬 Add Your Comment:"),
			dbc.Textarea(
				id="new-comment",
				placeholder="Enter your comment...",
				style={"minHeight": "100px"}
			),
			html.Div([
				html.Label("👤 Your Name:", className="mt-3 mb-1"),
				dbc.Input(
					id="user-name",
					value="Dr. Current",
					placeholder="Enter your name"
				)
			])
		]),
		dbc.ModalFooter([
			dbc.Button("Cancel", id="close-modal", className="me-2", color="secondary"),
			dbc.Button("💬 Add Comment", id="add-comment", color="primary")
		])
	], id="comment-modal", is_open=False, size="lg")

def create_loading_component():
	"""Create a loading spinner component"""
	return html.Div([
		dbc.Spinner([
			html.Div([
				html.I(className="fas fa-dna fa-3x text-primary mb-3"),
				html.H5("Loading variants...", className="text-primary")
			], className="text-center")
		], size="lg", color="primary")
	], className="loading-spinner")

def create_error_component(error_message):
	"""Create an error display component"""
	return html.Div([
		dbc.Alert([
			html.Div([
				html.I(className="fas fa-exclamation-triangle fa-2x text-danger mb-3"),
				html.H5("Error", className="text-danger mb-3"),
				html.P(error_message, className="mb-0")
			], className="text-center")
		], color="danger", className="border-0")
	])

def create_success_component(message):
	"""Create a success display component"""
	return html.Div([
		dbc.Alert([
			html.Div([
				html.I(className="fas fa-check-circle fa-2x text-success mb-3"),
				html.H5("Success", className="text-success mb-3"),
				html.P(message, className="mb-0")
			], className="text-center")
		], color="success", className="border-0")
	])

def create_no_selection_display():
	"""Create display for when no samples are selected"""
	return html.Div([
		dbc.Alert([
			html.Div([
				html.I(className="fas fa-users fa-2x text-warning mb-3"),
				html.H5("Please Select Samples", className="text-warning mb-3"),
				html.P("Select one or more samples to display their variants.", className="mb-0")
			], className="text-center")
		], color="warning", className="border-0", style={
			"background": "rgba(255, 193, 7, 0.1)", 
			"borderRadius": "15px", 
			"padding": "30px"
		})
	])

def create_variant_count_display(count, total_count, sample_count):
	"""Create variant count display"""
	sample_text = f" from {sample_count} sample{'s' if sample_count > 1 else ''}" if sample_count > 0 else ""
	
	return html.Div([
		html.H5([
			html.Span(f"{count:,}", className="text-primary fw-bold", style={"fontSize": "1.6rem"}),
			html.Span(" variants", style={"fontSize": "1.2rem"}),
			html.Span(sample_text, style={"fontSize": "1rem", "color": "#0097A7"})
		], className="mb-0"),
		html.Small(f"(from {total_count:,} total)", className="text-muted", 
				  style={"fontSize": "14px"}) if total_count != count else html.Div()
	])