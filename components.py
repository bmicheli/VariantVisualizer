"""
UI Components for Variant Visualizer
Contains all display components and layout functions
"""

import dash_bootstrap_components as dbc
from dash import html, dcc
import pandas as pd
import polars as pl
from database import get_variant_comments, db
from utils import *
from config import *

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

def create_beautiful_variant_display(df):
	"""Create optimized variant table display"""
	# Check DataFrame length properly for both Polars and Pandas
	if isinstance(df, pl.DataFrame):
		is_empty = len(df) == 0
	else:
		is_empty = df.empty
	
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
	
	# Limit display for performance (show first 1000 rows)
	display_df = df_pandas.head(MAX_DISPLAY_VARIANTS)
	
	# Create table rows efficiently
	table_rows = []
	
	for idx, variant in display_df.iterrows():
		variant_key = variant.get('variant_key', f"{variant['CHROM']}:{variant['POS']}:{variant['REF']}:{variant['ALT']}")
		
		# Format data
		position = f"{variant['CHROM']}:{variant['POS']:}"
		variant_text = f"{variant['REF']}>{variant['ALT']}"
		vaf_display = format_percentage(variant.get('VAF', 0))
		
		# Get population frequency (try multiple columns)
		pop_af = variant.get('af', variant.get('gnomad_af', variant.get('population_af', 0)))
		pop_af_display = format_frequency(pop_af)
		
		# Create data row
		data_row = html.Tr([
			# Sample ID
			html.Td(variant['SAMPLE'], style={"fontSize": "11px", "fontWeight": "500"}),
			# Position
			html.Td(position, style={"fontFamily": "monospace", "fontSize": "11px"}),
			# Gene
			html.Td(variant.get('gene', 'UNKNOWN'), style={"fontWeight": "bold", "color": "#0097A7", "fontSize": "12px"}),
			# Variant
			html.Td(variant_text, style={"fontFamily": "monospace", "fontSize": "11px", "backgroundColor": "#f8f9fa", "borderRadius": "4px", "textAlign": "center"}),
			# Genotype
			html.Td([get_genotype_badge(variant.get('GT', './.'))]),
			# VAF
			html.Td(vaf_display, style={"fontFamily": "monospace", "fontSize": "10px", "textAlign": "right"}),
			# Consequence
			html.Td([get_consequence_badge(variant.get('consequence', 'variant'))]),
			# ClinVar Classification
			html.Td([get_clinvar_badge(variant.get('clinvar_sig'))]),
			# Population Frequency (nouvelle colonne)
			html.Td([
				html.Span(pop_af_display, 
					style={
						"fontFamily": "monospace", 
						"fontSize": "10px", 
						"fontWeight": "bold",
						"color": "#dc3545" if pop_af and pop_af < 0.001 
							else "#ffc107" if pop_af and 0.001 <= pop_af < 0.01 
							else "#28a745" if pop_af and pop_af >= 0.01 
							else "#6c757d"
					},
					title=f"Population Allele Frequency: {pop_af_display}"
				)
			], style={"textAlign": "right"}),
			# Comments
			html.Td([
				dbc.Button([html.I(className="fas fa-comments me-1"), str(variant.get('comment_count', 0))], 
						id={"type": "comment-btn", "variant": variant_key, "sample": variant['SAMPLE']},
						color="outline-primary", size="sm", n_clicks=0)
			], style={"textAlign": "center"}),
		], className="variant-row", id={"type": "variant-row", "variant": variant_key, "sample": variant['SAMPLE']}, n_clicks=0)
		
		# Accordion collapse row for detailed view
		accordion_row = html.Tr([
			html.Td([
				dbc.Collapse([
					html.Div([
						create_variant_details_accordion(variant)
					], className="variant-details-content")
				], 
				id={"type": "variant-collapse", "variant": variant_key, "sample": variant['SAMPLE']},
				is_open=False)
			], colSpan=10, style={"padding": "0", "border": "none"})
		], style={"border": "none"})
		
		table_rows.extend([data_row, accordion_row])
	
	# Add notice if showing limited results
	total_count = len(df_pandas)
	showing_count = len(display_df)
	
	display_info = []
	if total_count > showing_count:
		display_info = [
			dbc.Alert([
				html.I(className="fas fa-info-circle me-2"),
				f"Showing first {showing_count:,} of {total_count:,} variants. Use filters to narrow results."
			], color="info", className="mb-3")
		]
	
	# Create table with headers
	return html.Div([
		*display_info,
		html.Table([
			html.Thead([
				html.Tr([
					html.Th("Sample", className="sortable-header"),
					html.Th("Position", className="sortable-header"),
					html.Th("Gene", className="sortable-header"),
					html.Th("Variant", className="sortable-header"),
					html.Th("Genotype", className="sortable-header"),
					html.Th("VAF", className="sortable-header"),
					html.Th("Consequence", className="sortable-header"),
					html.Th("ClinVar", className="sortable-header"),
					html.Th("Pop. AF", className="sortable-header", title="Population Allele Frequency"),
					html.Th("Comments", className="sortable-header"),
				])
			]),
			html.Tbody(table_rows)
		], className="variants-table")
	], className="variants-table-container")

def create_variant_details_accordion(variant):
	"""Create detailed information panel for accordion expansion"""
	
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
							html.Span(f"{variant['CHROM']}:{variant['POS']:}")
						], className="mb-2"),
						
						html.Div([
							html.Strong("Reference: "),
							html.Span(variant['REF'])
						], className="mb-2"),
						
						html.Div([
							html.Strong("Alternate: "),
							html.Span(variant['ALT'])
						], className="mb-2"),
						
						html.Div([
							html.Strong("Genotype: "),
							get_genotype_badge(variant.get('GT', './.'))
						], className="mb-2"),
						
						html.Div([
							html.Strong("Quality: "),
							html.Span(f"{variant.get('QUAL', 0):.1f}" if pd.notna(variant.get('QUAL')) else "N/A")
						], className="mb-2"),
						
						# Ajout du statut de review dans les détails
						html.Div([
							html.Strong("Review Status: "),
							get_status_badge(variant.get('review_status', 'Pending'))
						], className="mb-2")
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
							html.Span(format_percentage(variant.get('VAF', 0)))
						], className="mb-2"),
						
						html.Div([
							html.Strong("Depth: "),
							html.Span(f"{variant.get('DP', 0)}" if pd.notna(variant.get('DP')) else "N/A")
						], className="mb-2"),
						
						html.Div([
							html.Strong("Allelic Depth: "),
							html.Span(variant.get('AD', 'N/A'))
						], className="mb-2"),
						
						html.Div([
							html.Strong("Genotype Quality: "),
							html.Span(f"{variant.get('GQ', 0):.1f}" if pd.notna(variant.get('GQ')) else "N/A")
						], className="mb-2"),
						
						html.Div([
							html.Strong("Sample ID: "),
							html.Span(sample_id)
						], className="mb-2")
					])
				], className="detail-section uniform-height")
			], width=4),
			
			# Column 3: Clinical Information
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
						
						html.Div([
							html.Strong("Gene: "),
							html.Span(variant.get('gene', 'UNKNOWN'))
						], className="mb-2"),
						
						html.Div([
							html.Strong("Consequence: "),
							get_consequence_badge(variant.get('consequence', 'variant'))
						], className="mb-2"),
						
						html.Div([
							html.Strong("AA Change: "),
							html.Span(variant.get('aa_change', 'N/A'))
						], className="mb-2"),
						
						html.Div([
							html.Strong("gnomAD AF: "),
							html.Span(format_frequency(variant.get('gnomad_af', 0)))
						], className="mb-2")
					])
				], className="detail-section uniform-height")
			], width=4)
		]),
		
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
									else "#6c757d"}
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
									else "#6c757d"}
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
									else "#6c757d"}
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
									else "#6c757d"}
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
									else "#6c757d"}
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
									else "#6c757d"}
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
						"Population & Other Scores"
					], className="text-primary mb-3"),
					
					html.Div([
						html.Div([
							html.Strong("PrimateAI Score: "),
							html.Span(
								format_score(variant.get('primateai_score', 0)) if variant.get('primateai_score') else "N/A",
								className="fw-bold",
								style={"color": "#28a745" if variant.get('primateai_score', 0) and variant.get('primateai_score', 0) < 0.5 
									else "#ffc107" if variant.get('primateai_score', 0) and 0.5 <= variant.get('primateai_score', 0) < 0.8 
									else "#dc3545" if variant.get('primateai_score', 0) and variant.get('primateai_score', 0) >= 0.8 
									else "#6c757d"}
							),
							html.Br(),
							html.Small(" (>0.8 pathogenic)", className="text-muted")
						], className="mb-2"),
						
						html.Div([
							html.Strong("Population AF: "),
							html.Span(
								format_frequency(variant.get('af', 0)) if variant.get('af') else "N/A",
								className="fw-bold",
								style={"color": "#dc3545" if variant.get('af', 0) and variant.get('af', 0) < 0.001 
									else "#ffc107" if variant.get('af', 0) and 0.001 <= variant.get('af', 0) < 0.01 
									else "#28a745" if variant.get('af', 0) and variant.get('af', 0) >= 0.01 
									else "#6c757d"}
							),
							html.Br(),
							html.Small(" (rare <0.1%)", className="text-muted")
						], className="mb-2"),
						
						html.Div([
							html.Strong("Allele Count: "),
							html.Span(f"{variant.get('ac', 0)}" if pd.notna(variant.get('ac')) else "N/A"),
							html.Br(),
							html.Small(" (rare <0.1%)", className="text-muted")
						], className="mb-2")
					])
				], className="detail-section uniform-height")
			], width=4)
		]) if any([variant.get('cadd_score'), variant.get('sift_score'), variant.get('polyphen_score'), 
				variant.get('revel_score'), variant.get('splice_ai'), variant.get('pli_score'), 
				variant.get('primateai_score'), variant.get('af')]) else html.Div(),
			
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
		])
	], style={"padding": "20px", "backgroundColor": "rgba(240, 248, 255, 0.8)", "borderRadius": "0 0 12px 12px"})

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
	"""Create the advanced filters sidebar with corrected default values"""
	return html.Div([
		# Sidebar header
		html.Div([
			html.H4([
				html.I(className="fas fa-sliders-h me-2"),
				"Advanced Filters"
			], className="text-primary mb-0"),
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
				html.Label("🔍 Search", className="fw-bold mb-2"),
				dbc.Input(
					id="search-input",
					placeholder="Search genes, samples...",
					debounce=True
				)
			], className="mb-4"),
			
			# Custom Filters
			html.Div([
				html.Label("⚙️ Advanced Filters", className="fw-bold mb-3"),
				
				html.Label("Genotype:", className="small mb-1"),
				dcc.Dropdown(
					id="genotype-filter",
					options=GENOTYPE_OPTIONS,
					value="all",
					clearable=False,
					className="mb-3"
				),
				
				html.Label("Chromosome:", className="small mb-1"),
				dcc.Dropdown(
					id="chromosome-filter",
					options=CHROMOSOME_OPTIONS,
					value="all",
					clearable=False,
					className="mb-3"
				),
				
				html.Label("VAF Range:", className="small mb-1"),
				dcc.RangeSlider(
					id="vaf-range-slider",
					min=0,
					max=1,
					step=0.01,
					value=[0, 1],
					marks={0: '0%', 0.25: '25%', 0.5: '50%', 0.75: '75%', 1: '100%'},
					tooltip={"placement": "bottom", "always_visible": True}
				),
				
				html.Label("Quality Score (min):", className="small mb-1 mt-3"),
				dbc.Input(
					id="quality-filter",
					type="number",
					placeholder="Minimum quality score",
					min=0,
					step=0.1,
					className="mb-3"
				),
				
				# CORRIGÉ: Valeurs par défaut plus permissives pour Population AF
				html.Label("Population AF Range:", className="small mb-1 mt-3"),
				dcc.RangeSlider(
					id="pop-af-range-slider",
					min=0,
					max=1,  # Changé de 0.5 à 1.0 (100%)
					step=0.001,
					value=[0, 1],  # Changé de [0, 0.5] à [0, 1]
					marks={
						0: '0%', 
						0.001: '0.1%', 
						0.01: '1%', 
						0.1: '10%', 
						0.5: '50%',
						1: '100%'  # Ajouté
					},
					tooltip={"placement": "bottom", "always_visible": True}
				),
				
				# Note d'aide pour l'utilisateur
				html.Small([
					html.I(className="fas fa-info-circle me-1"),
					"Population AF: Leave at 0-100% to include all variants. Lower values = rarer variants."
				], className="text-muted mt-1 d-block")
			]),
			
			# Clear Filters
			dbc.Button(
				"🗑️ Clear All Filters", 
				id="clear-filters", 
				color="outline-danger", 
				size="sm", 
				className="w-100 mt-3"
			)
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
					], className="mb-0"),
				], style={"flex": "1"}),
				html.Div([
					dbc.ButtonGroup([
						dbc.Button([
							html.I(className="fas fa-sync-alt me-2"),
							"Refresh Data"
						], id="refresh-data-btn", color="info", outline=True, className="me-2"),
						dbc.Button([
							html.I(className="fas fa-download me-2"),
							"Export"
						], id="export-btn", color="success", outline=True),
						dbc.Button([
							html.I(className="fas fa-chart-bar me-2"),
							"Stats"
						], id="stats-btn", color="warning", outline=True, className="ms-2"),
					])
				])
			], style={"display": "flex", "alignItems": "center", "justifyContent": "space-between"})
		])
	], className="glass-card mb-3")

def create_sample_selector():
	"""Create the sample selector component - FIXED Z-INDEX"""
	return dbc.Card([
		dbc.CardBody([
			html.Div([
				html.Div([
					html.H6([html.I(className="fas fa-users me-2"), "Sample Selection"], className="mb-2 text-primary"),
					html.P("Select samples to display variants:", className="text-muted mb-0 small")
				]),
				html.Div([
					dcc.Dropdown(
						id="sample-selector",
						placeholder="Select samples...",
						multi=True,
						searchable=True,
						clearable=True,
						style={"minWidth": "300px"}
					)
				], style={"flex": "1", "position": "relative"}),
				dbc.ButtonGroup([
					dbc.Button("Clear", id="clear-samples", color="outline-secondary", size="sm")
				])
			], style={"display": "flex", "alignItems": "center", "gap": "15px"})
		], style={"padding": "20px"})
	], className="sample-selector-container mb-3")

def create_main_filters_panel():
	"""Create the main filters panel with quick filters"""
	return dbc.Card([
		dbc.CardBody([
			dbc.Row([
				dbc.Col([
					html.Div(id="variant-count")
				], width=3),
				dbc.Col([
					html.Div([
						html.Label("⚡ Quick Filters:", className="fw-bold me-3 mb-0"),
						html.Div([
							dbc.Button([f["icon"], " ", f["name"]], 
									id={"type": "preset-filter", "id": f["id"]},
									color="outline-info", size="sm", className="quick-filter-btn me-2 mb-1", n_clicks=0)
							for f in PRESET_FILTERS
						], style={"display": "flex", "flexWrap": "wrap"})
					])
				], width=7),
				dbc.Col([
					dbc.Button("More Filters", id="more-filters-btn", color="outline-secondary", size="sm")
				], width=2)
			], align="center")
		], style={"padding": "15px"})
	], className="main-filters-panel mb-3")

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

def create_stats_modal():
	"""Create statistics modal component"""
	return dbc.Modal([
		dbc.ModalHeader([
			dbc.ModalTitle([
				html.I(className="fas fa-chart-bar me-2"),
				"Database Statistics"
			])
		]),
		dbc.ModalBody([
			html.Div(id="stats-content")
		]),
		dbc.ModalFooter([
			dbc.Button("Close", id="close-stats-modal", color="secondary")
		])
	], id="stats-modal", is_open=False, size="xl")

def create_stats_display(stats):
	"""Create statistics display content"""
	if not stats or stats.get('total_variants', 0) == 0:
		return html.Div([
			dbc.Alert("No statistics available - database is empty.", color="info")
		])
	
	return html.Div([
		dbc.Row([
			# Overall Statistics
			dbc.Col([
				dbc.Card([
					dbc.CardHeader("📊 Overall Statistics"),
					dbc.CardBody([
						html.P([html.Strong("Total Variants: "), f"{stats.get('total_variants', 0):,}"]),
						html.P([html.Strong("Unique Genes: "), f"{stats.get('unique_genes', 0):,}"]),
						html.P([html.Strong("Samples: "), f"{stats.get('samples_with_variants', 0):,}"]),
						html.P([html.Strong("Reviewed: "), f"{stats.get('reviewed_variants', 0):,}"]),
						html.P([html.Strong("ClinVar Annotated: "), f"{stats.get('clinvar_annotated', 0):,}"])
					])
				])
			], width=6),
			
			# Distribution Charts would go here
			dbc.Col([
				dbc.Card([
					dbc.CardHeader("🧬 Consequence Distribution"),
					dbc.CardBody([
						html.Div(id="consequence-chart")
					])
				])
			], width=6)
		], className="mb-3"),
		
		dbc.Row([
			dbc.Col([
				dbc.Card([
					dbc.CardHeader("🧬 Chromosome Distribution"),
					dbc.CardBody([
						html.Div(id="chromosome-chart")
					])
				])
			], width=6),
			
			dbc.Col([
				dbc.Card([
					dbc.CardHeader("📋 ClinVar Distribution"),
					dbc.CardBody([
						html.Div(id="clinvar-chart")
					])
				])
			], width=6)
		])
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
			html.Span(f"{count:,}", className="text-primary fw-bold", style={"fontSize": "1.5em"}),
			html.Span(" variants", style={"fontSize": "1.2em"}),
			html.Span(sample_text, style={"fontSize": "0.9em", "color": "#0097A7"})
		], className="mb-0"),
		html.Small(f"(from {total_count:,} total)", className="text-muted") if total_count != count else html.Div()
	])