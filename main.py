import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, Output, Input, State, callback_context, ALL, MATCH
import pandas as pd
import json
from datetime import datetime
# Import local modules
from config import *
from database import *
from utils import *
from components import *
from gene_panels import (init_gene_panels, update_panels_if_needed, get_available_panels,get_genes_for_panels, get_panel_info, force_update_panels)

# =============================================================================
# APP INITIALIZATION
# =============================================================================

# Initialize Parquet database
print("🔧 Initializing Parquet database...")
init_parquet_database()
print("✅ Parquet database structure initialized.")

# Initialize gene panels
print("🧬 Initializing gene panel system...")
init_gene_panels()
update_panels_if_needed()
print("✅ Gene panel system initialized.")

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=EXTERNAL_STYLESHEETS, suppress_callback_exceptions=True)

# Set app title
app.title = "Variant Visualizer"

# Add custom CSS
app.index_string = f'''
<!DOCTYPE html>
<html>
    <head>
        {{%metas%}}
        <title>{{%title%}}</title>
        {{%favicon%}}
        {{%css%}}
        <style>
            {CUSTOM_CSS}
        </style>
    </head>
    <body>
        {{%app_entry%}}
        <footer>
            {{%config%}}
            {{%scripts%}}
            {{%renderer%}}
        </footer>
    </body>
</html>
'''

# =============================================================================
# APP LAYOUT
# =============================================================================

app.layout = html.Div([
    # Sidebar overlay
    html.Div(id="sidebar-overlay", className="sidebar-overlay"),
    
    # Sidebar with advanced filters
    create_sidebar(),
    
    # Main container
    dbc.Container([
        # Header
        create_header(),
        
        # Sample Selector Panel
        create_sample_selector(),
        
        # Gene Panel Selector Panel
        create_gene_panel_selector(),
        
        # Main Filters Panel (with reset button)
        create_main_filters_panel(),
        
        # Variants Display
        dbc.Card([
            dbc.CardBody([
                html.Div(id="variants-display")
            ])
        ], className="glass-card")
        
    ], fluid=True, style={"padding": "20px"}),
    
    # Modals
    create_comment_modal(),
    create_panel_info_modal(),
    
    # Toast notifications
    create_update_status_toast(),
    
    # Data stores
    dcc.Store(id="filtered-variants"),
    dcc.Store(id="selected-variant-id"),
    dcc.Store(id="active-filters", data={}),
    dcc.Store(id="sidebar-open", data=False),
    dcc.Store(id="available-samples", data=[]),
    dcc.Store(id="selected-gene-panels", data=[]),
    dcc.Store(id="panel-genes", data=[]),
    dcc.Store(id="sort-state", data={"column": None, "direction": "asc"})  # NOUVEAU: Store pour le tri
    
], style={"minHeight": "100vh", "background": "linear-gradient(135deg, #00BCD4 0%, #4DD0E1 50%, #80E5A3 100%)"})

# =============================================================================
# CALLBACKS
# =============================================================================

# Sample options callback
@app.callback(
    [Output("sample-selector", "options"), Output("available-samples", "data")],
    [Input("sample-selector", "id")]
)
def update_sample_options(sample_selector_id):
    """Update sample selector options when app loads"""
    samples = get_available_samples()
    options = [{"label": sample, "value": sample} for sample in samples]
    return options, samples

# Gene panel options callback
@app.callback(
    Output("gene-panel-selector", "options"),
    [Input("gene-panel-selector", "id"), Input("update-panels-btn", "n_clicks")]
)
def update_gene_panel_options(panel_selector_id, update_clicks):
    """Update gene panel selector options"""
    return get_available_panels()

# Gene panel selection callback 
@app.callback(
    [Output("selected-gene-panels", "data"), 
     Output("panel-genes", "data")],
    [Input("gene-panel-selector", "value"), 
     Input("clear-gene-panels", "n_clicks")],
    [State("selected-gene-panels", "data")],
    prevent_initial_call=True
)
def handle_gene_panel_selection(selected_panels, clear_clicks, current_panels):
    """Handle gene panel selection - SIMPLIFIED (no green gene filtering)"""
    ctx = callback_context
    if not ctx.triggered:
        return current_panels or [], []
    
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    
    # Clear panels
    if trigger == "clear-gene-panels":
        return [], []
    
    # Get genes for selected panels
    if selected_panels:
        genes = get_genes_for_panels(selected_panels)
        logger.info(f"Selected panels: {selected_panels}, Genes: {len(genes)}")
    else:
        genes = []
    
    return selected_panels or [], genes

# Callback pour la gestion du tri des colonnes
@app.callback(
    Output("sort-state", "data"),
    [Input({"type": "sort-header", "column": ALL}, "n_clicks")],
    [State("sort-state", "data")],
    prevent_initial_call=True
)
def handle_column_sort(n_clicks_list, current_sort):
    """Handle column sorting when headers are clicked"""
    ctx = callback_context
    if not ctx.triggered or not any(n_clicks_list):
        return current_sort
    
    # Get the triggered button ID
    trigger_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])
    clicked_column = trigger_id['column']
    
    # Determine new sort direction
    if current_sort['column'] == clicked_column:
        # Same column clicked - toggle direction
        new_direction = "desc" if current_sort['direction'] == "asc" else "asc"
    else:
        # Different column clicked - default to ascending
        new_direction = "asc"
    
    return {"column": clicked_column, "direction": new_direction}

# Panel info modal callback
@app.callback(
    [Output("panel-info-modal", "is_open"), Output("panel-info-content", "children")],
    [Input("panel-info-btn", "n_clicks"), Input("close-panel-info-modal", "n_clicks")],
    [State("panel-info-modal", "is_open"), State("selected-gene-panels", "data")]
)
def handle_panel_info_modal(info_clicks, close_clicks, is_open, selected_panels):
    """Handle panel info modal"""
    ctx = callback_context
    if not ctx.triggered:
        return False, html.Div()
    
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger == "panel-info-btn" and selected_panels:
        content = create_panel_info_content(selected_panels, selected_panels)
        return True, content
    elif trigger == "close-panel-info-modal":
        return False, html.Div()
    
    return is_open, html.Div()

# Update panels callback
@app.callback(
    [Output("update-status-toast", "is_open"), Output("update-status-toast", "children")],
    [Input("update-panels-btn", "n_clicks")],
    prevent_initial_call=True
)
def handle_panel_update(update_clicks):
    """Handle panel update button"""
    if update_clicks:
        try:
            force_update_panels()
            message = html.Div([
                html.I(className="fas fa-check-circle me-2"),
                "Gene panels updated successfully!"
            ])
            return True, message
        except Exception as e:
            message = html.Div([
                html.I(className="fas fa-exclamation-triangle me-2"),
                f"Error updating panels: {str(e)}"
            ])
            return True, message
    
    return False, html.Div()

# Sample selection callback
@app.callback(
    Output("sample-selector", "value"),
    [Input("clear-samples", "n_clicks")],
    prevent_initial_call=True
)
def handle_clear_samples(clear_clicks):
    """Handle clear samples button"""
    if clear_clicks:
        return []
    return dash.no_update

# Gene panel selection callback
@app.callback(
    Output("gene-panel-selector", "value"),
    [Input("clear-gene-panels", "n_clicks")],
    prevent_initial_call=True
)
def handle_clear_gene_panels(clear_clicks):
    """Handle clear gene panels button"""
    if clear_clicks:
        return []
    return dash.no_update

# Updated callback for sidebar with apply button functionality
@app.callback(
    [Output("filter-sidebar", "className"), Output("sidebar-overlay", "className"), Output("sidebar-open", "data")],
    [Input("more-filters-btn", "n_clicks"), 
     Input("close-sidebar-btn", "n_clicks"), 
     Input("sidebar-overlay", "n_clicks"),
     Input("apply-filters-btn", "n_clicks")],
    [State("sidebar-open", "data")],
    prevent_initial_call=True
)
def toggle_sidebar(toggle_clicks, close_clicks, overlay_clicks, apply_clicks, is_open):
    """Toggle sidebar visibility - closes when apply is clicked"""
    ctx = callback_context
    if not ctx.triggered:
        return "filter-sidebar", "sidebar-overlay", False
    
    trigger = ctx.triggered[0]['prop_id'].split('.')[0]
    
    # Close sidebar when apply button is clicked
    if trigger == "apply-filters-btn":
        return "filter-sidebar", "sidebar-overlay", False
    
    # Toggle logic for other buttons
    new_state = not is_open
    sidebar_class = "filter-sidebar open" if new_state else "filter-sidebar"
    overlay_class = "sidebar-overlay open" if new_state else "sidebar-overlay"
    return sidebar_class, overlay_class, new_state

# Callback principal avec tri intégré - FIXED GENE SORTING
@app.callback(
    [Output("variants-display", "children"), 
     Output("variant-count", "children"), 
     Output("filtered-variants", "data")],
    [Input("apply-filters-btn", "n_clicks"),  
     Input("sample-selector", "value"),       
     Input("active-filters", "data"),         
     Input("panel-genes", "data"),
     Input("sort-state", "data")],
    [State("search-input", "value"), 
     State("vaf-range-slider", "value"),
     State("selected-gene-panels", "data")],     
    prevent_initial_call=False
)
def update_variants_display_optimized(apply_clicks, selected_samples, active_filters, panel_genes, 
                                      sort_state, search_term, vaf_range, selected_panels):
    """SUPER OPTIMIZED variants display with HGNC_ID-aware filtering and SORTING - FIXED GENE SORTING"""
    
    if not selected_samples:
        no_selection = create_no_selection_display()
        count_display = create_variant_count_display(0, 0, 0)
        return no_selection, count_display, []
    
    try:
        import time
        start_time = time.time()
        
        # Charger la DB
        df = db.load_variants_lazy(samples=selected_samples, limit=MAX_LOAD_LIMIT)
        if is_dataframe_empty(df):
            empty_display = create_beautiful_variant_display(df)
            empty_count = create_variant_count_display(0, 0, len(selected_samples))
            return empty_display, empty_count, []
        
        original_count = len(df)
        
        # === FILTRE PANEL DE GENES ===
        if panel_genes:
            logger.info(f"Applying gene panel filter with {len(panel_genes)} entries")

            # Charger mapping HGNC_ID → symbol
            gene_mapping = load_gene_mapping() or {}
            reverse_mapping = {v.upper(): k for k, v in gene_mapping.items()}  # symbol → ID

            # Construire l'ensemble des HGNC_ID à filtrer
            panel_ids = set()
            for gene in panel_genes:
                g_upper = str(gene).upper()
                if g_upper in reverse_mapping:
                    panel_ids.add(reverse_mapping[g_upper])   # symbole → ID
                else:
                    panel_ids.add(g_upper)  # déjà HGNC_ID ou inconnu
            
            logger.info(f"Filtering variants on {len(panel_ids)} HGNC_IDs")

            # Conversion en Polars si besoin
            if isinstance(df, pd.DataFrame):
                df = pl.from_pandas(df)

            gene_list = list(panel_ids)

            # Si plusieurs IDs dans la colonne "gene" (ex: "8654,8654"), on split et on check la présence
            if len(gene_list) == 1:
                df = df.filter(pl.col("gene").str.split(",").list.contains(gene_list[0]))
            else:
                conditions = [pl.col("gene").str.split(",").list.contains(g) for g in gene_list]
                df = df.filter(pl.any_horizontal(conditions))

            logger.info(f"After gene panel filter: {len(df)} variants (from {original_count})")
        
        # === AUTRES FILTRES ===
        df = apply_filters_optimized(
            df, 
            search_term=search_term,
            genotype_filter=None,
            chromosome_filter=None,
            active_filters=active_filters,
            selected_samples=selected_samples
        )
        
        # Filtre VAF
        if vaf_range and len(vaf_range) == 2:
            before_vaf = len(df)
            df = df.filter((pl.col('VAF') >= vaf_range[0]) & (pl.col('VAF') <= vaf_range[1]))
            logger.info(f"VAF filter {vaf_range}: {before_vaf} -> {len(df)}")
        
        # === TRI DES DONNÉES - FIXED GENE SORTING ===
        if sort_state and sort_state.get('column'):
            sort_column = sort_state['column']
            sort_direction = sort_state['direction']
            
            # Conversion en Polars si nécessaire
            if isinstance(df, pd.DataFrame):
                df = pl.from_pandas(df)
            
            # Appliquer le tri selon la colonne
            ascending = sort_direction == "asc"
            
            try:
                if sort_column == "sample":
                    df = df.sort("SAMPLE", descending=not ascending)
                elif sort_column == "position":
                    # Tri spécial pour position (chromosome puis position)
                    df = df.with_columns([
                        pl.when(pl.col('CHROM') == 'X').then(23)
                        .when(pl.col('CHROM') == 'Y').then(24)
                        .when(pl.col('CHROM') == 'MT').then(25)
                        .otherwise(pl.col('CHROM').cast(pl.Int32, strict=False))
                        .alias('chrom_numeric')
                    ])
                    df = df.sort(['chrom_numeric', 'POS'], descending=[not ascending, not ascending])
                    df = df.drop('chrom_numeric')
                elif sort_column == "gene":
                    # FIXED: Tri alphabétique par symbole de gène au lieu du HGNC ID
                    logger.info("Applying alphabetical gene symbol sorting")
                    
                    # Charger le mapping HGNC ID → symbole
                    gene_mapping = load_gene_mapping() or {}
                    
                    # Créer une colonne temporaire avec les symboles de gènes pour le tri
                    def convert_gene_ids_to_symbols(gene_str):
                        """Convert gene IDs to symbols for sorting"""
                        if pd.isna(gene_str) or gene_str in [None, '', 'UNKNOWN']:
                            return 'UNKNOWN'
                        
                        # Handle multiple genes separated by commas
                        genes = [g.strip() for g in str(gene_str).split(',') if g.strip()]
                        if not genes:
                            return 'UNKNOWN'
                        
                        # Convert each gene ID to symbol and take the first one for sorting
                        symbols = []
                        for gene in genes:
                            if gene in gene_mapping:
                                symbols.append(gene_mapping[gene])
                            else:
                                symbols.append(gene)  # Keep original if not found
                        
                        # Return first symbol for sorting (alphabetically)
                        return symbols[0] if symbols else 'UNKNOWN'
                    
                    # Apply the conversion to create a temporary sorting column
                    df = df.with_columns([
                        pl.col('gene').map_elements(convert_gene_ids_to_symbols, return_dtype=pl.Utf8).alias('gene_symbol_for_sort')
                    ])
                    
                    # Sort by the gene symbol column
                    df = df.sort("gene_symbol_for_sort", descending=not ascending)
                    
                    # Remove the temporary column
                    df = df.drop('gene_symbol_for_sort')
                    
                    logger.info(f"Applied alphabetical gene sorting: {sort_direction}")
                elif sort_column == "genotype":
                    df = df.sort("GT", descending=not ascending)
                elif sort_column == "moi":
                    df = df.sort("moi", descending=not ascending)
                elif sort_column == "vaf":
                    df = df.sort("VAF", descending=not ascending)
                elif sort_column == "consequence":
                    df = df.sort("consequence", descending=not ascending)
                elif sort_column == "clinvar":
                    df = df.sort("clinvar_sig", descending=not ascending)
                elif sort_column == "gnomad_af":
                    # Tri par max_gnomad_af si disponible, sinon gnomad_af
                    if 'max_gnomad_af' in df.columns:
                        df = df.sort("max_gnomad_af", descending=not ascending)
                    else:
                        df = df.sort("gnomad_af", descending=not ascending)
                
                logger.info(f"Applied sort: {sort_column} {sort_direction}")
            except Exception as e:
                logger.error(f"Error applying sort: {e}")
        
        # === AFFICHAGE & STOCKAGE ===
        display = create_beautiful_variant_display(
            df, 
            sort_column=sort_state.get('column') if sort_state else None,
            sort_direction=sort_state.get('direction') if sort_state else "asc"
        )
        count_display = create_variant_count_display(len(df), original_count, len(selected_samples))
        
        try:
            if isinstance(df, pl.DataFrame):
                storage_data = df.head(MAX_DISPLAY_VARIANTS).to_pandas().to_dict('records')
            else:
                storage_data = df.head(MAX_DISPLAY_VARIANTS).to_dict('records')
        except Exception as e:
            logger.error(f"Error converting to storage format: {e}")
            storage_data = []
        
        return display, count_display, storage_data
    
    except Exception as e:
        logger.error(f"Error in update_variants_display: {e}")
        error_display = create_error_component(f"Error loading variants: {str(e)}")
        return error_display, create_variant_count_display(0, 0, 0), []

# NEW: Lazy loading callback for variant details
@app.callback(
    Output({"type": "variant-details-lazy", "variant": MATCH, "sample": MATCH}, "children"),
    [Input({"type": "variant-collapse", "variant": MATCH, "sample": MATCH}, "is_open")],
    [State("filtered-variants", "data")],
    prevent_initial_call=True
)
def load_variant_details_lazy(is_open, filtered_variants):
    """Lazy load variant details only when accordion is opened"""
    if not is_open or not filtered_variants:
        return html.Div()
    
    ctx = callback_context
    if not ctx.triggered:
        return html.Div()
    
    # Extract variant and sample from triggered component
    trigger_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])
    variant_key = trigger_id.get('variant')
    sample_id = trigger_id.get('sample')
    
    # Find the variant in filtered data
    for variant_dict in filtered_variants:
        if (variant_dict.get('variant_key') == variant_key and 
            variant_dict.get('SAMPLE') == sample_id):
            # Ensure max_gnomad_af is available for display
            if 'max_gnomad_af' not in variant_dict or pd.isna(variant_dict.get('max_gnomad_af')):
                variant_dict['max_gnomad_af'] = get_max_gnomad_af_from_variant(variant_dict)
            
            return create_variant_details_accordion(pd.Series(variant_dict))
    
    return html.Div([html.P("Variant details not found.")])

# NEW: AA Change toggle callback
@app.callback(
    Output({"type": "aa-change-collapse", "variant": MATCH, "sample": MATCH}, "is_open"),
    [Input({"type": "aa-change-toggle", "variant": MATCH, "sample": MATCH}, "n_clicks")],
    [State({"type": "aa-change-collapse", "variant": MATCH, "sample": MATCH}, "is_open")],
    prevent_initial_call=True
)
def toggle_aa_change_details(n_clicks, is_open):
    """Toggle AA change details collapse"""
    if n_clicks:
        return not is_open
    return is_open

# Reset all filters callback - MAIN PANEL RESET BUTTON
@app.callback(
    [Output("search-input", "value", allow_duplicate=True), 
     Output("active-filters", "data", allow_duplicate=True),
     Output("vaf-range-slider", "value", allow_duplicate=True),
     Output("sort-state", "data", allow_duplicate=True)],  # NOUVEAU: Reset du tri
    [Input("reset-all-btn", "n_clicks")],
    prevent_initial_call=True
)
def reset_all_filters(n_clicks):
    """Reset all filters from the main panel reset button"""
    if n_clicks:
        return "", {}, [0, 1], {"column": None, "direction": "asc"}
    return dash.no_update

# Updated Clear filters callback to include apply functionality
@app.callback(
    [Output("search-input", "value", allow_duplicate=True), 
     Output("active-filters", "data", allow_duplicate=True),
     Output("vaf-range-slider", "value", allow_duplicate=True),
     Output("filter-sidebar", "className", allow_duplicate=True), 
     Output("sidebar-overlay", "className", allow_duplicate=True), 
     Output("sidebar-open", "data", allow_duplicate=True),
     Output("sort-state", "data", allow_duplicate=True)],  # NOUVEAU: Reset du tri
    [Input("clear-filters", "n_clicks")],
    prevent_initial_call=True
)
def clear_all_filters_sidebar(n_clicks):
    """Clear filters from sidebar clear button and close sidebar"""
    if n_clicks:
        return "", {}, [0, 1], "filter-sidebar", "sidebar-overlay", False, {"column": None, "direction": "asc"}
    return dash.no_update

# Variant row expansion callback - RETOUR À LA VERSION ORIGINALE
@app.callback(
    Output({"type": "variant-collapse", "variant": MATCH, "sample": MATCH}, "is_open"),
    [Input({"type": "variant-row", "variant": MATCH, "sample": MATCH}, "n_clicks")],
    [State({"type": "variant-collapse", "variant": MATCH, "sample": MATCH}, "is_open")]
)
def handle_variant_accordion(n_clicks, is_open):
    """Handle accordion expansion/collapse of variant rows"""
    if n_clicks:
        return not is_open
    return is_open

# Quick filters callback
@app.callback(
    Output("active-filters", "data"),
    [Input({"type": "preset-filter", "id": ALL}, "n_clicks")],
    [State("active-filters", "data")]
)
def handle_preset_filters(n_clicks_list, current_filters):
    """Handle preset filter buttons"""
    ctx = callback_context
    if not ctx.triggered or not any(n_clicks_list):
        return current_filters or {}
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    filter_data = json.loads(button_id)
    filter_id = filter_data['id']
    
    new_filters = current_filters.copy() if current_filters else {}
    
    if new_filters.get(filter_id):
        del new_filters[filter_id]
    else:
        new_filters[filter_id] = True
    
    return new_filters

# Comment modal callback
@app.callback(
    [Output("comment-modal", "is_open"), Output("comment-modal-title", "children"), 
     Output("existing-comments", "children"), Output("selected-variant-id", "data")],
    [Input({"type": "comment-btn", "variant": ALL, "sample": ALL}, "n_clicks"),
     Input({"type": "add-comment-btn", "variant": ALL, "sample": ALL}, "n_clicks"),
     Input("close-modal", "n_clicks"), Input("add-comment", "n_clicks")],
    [State("comment-modal", "is_open"), State("new-comment", "value"), 
     State("user-name", "value"), State("selected-variant-id", "data")]
)
def handle_comment_modal(comment_clicks, add_comment_clicks, close_clicks, add_clicks, 
                        is_open, new_comment, user_name, selected_variant_id):
    """Handle comment modal interactions"""
    ctx = callback_context
    if not ctx.triggered:
        return False, dash.no_update, dash.no_update, selected_variant_id

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

    try:
        button_data = json.loads(trigger_id)
        if button_data.get('type') in ['comment-btn', 'add-comment-btn'] and ctx.triggered[0]['value'] > 0:
            variant_key = button_data.get('variant')
            sample_id = button_data.get('sample')
            title = f"💬 Comments for {variant_key} - {sample_id}"
            comments_df = get_variant_comments(variant_key, sample_id)
            
            if not comments_df.empty:
                comments_display = [
                    dbc.Card([
                        dbc.CardBody([
                            html.H6(comment['user_name'], className="text-primary mb-1"),
                            html.P(comment['comment_text'], className="mb-1"),
                            html.Small(comment['timestamp'], className="text-muted")
                        ])
                    ], className="mb-2")
                    for _, comment in comments_df.iterrows()
                ]
            else:
                comments_display = [dbc.Alert("No comments yet.", color="info")]
            
            return True, title, comments_display, {"variant_key": variant_key, "sample_id": sample_id}
    except:
        pass

    if trigger_id == "add-comment" and new_comment and user_name and selected_variant_id:
        success = add_variant_comment(
            selected_variant_id["variant_key"], 
            selected_variant_id["sample_id"], 
            user_name, 
            new_comment
        )
        if success:
            return False, "", "", None
        else:
            return True, dash.no_update, dash.no_update, selected_variant_id

    if trigger_id == "close-modal":
        return False, "", "", None

    return is_open, dash.no_update, dash.no_update, selected_variant_id

# Clear comment field callback
@app.callback(
    Output("new-comment", "value"),
    [Input("add-comment", "n_clicks")],
    [State("new-comment", "value")]
)
def clear_comment_field(n_clicks, comment_value):
    """Clear comment field after adding comment"""
    if n_clicks and comment_value:
        return ""
    return dash.no_update

# SIMPLE: Callback clientside pour la sélection de texte uniquement
app.clientside_callback(
    """
    function() {
        setTimeout(() => {
            // Solution simple : empêcher seulement les événements de sélection de texte 
            // de déclencher les callbacks d'accordéon
            document.addEventListener('selectstart', function(e) {
                e.stopPropagation();
            });
            
            document.addEventListener('mouseup', function(e) {
                if (window.getSelection && window.getSelection().toString()) {
                    e.stopPropagation();
                }
            });
            
        }, 1000);
        
        return window.dash_clientside.no_update;
    }
    """,
    Output("variants-display", "style", allow_duplicate=True),
    [Input("variants-display", "children")],
    prevent_initial_call=True
)

# Force dropdown z-index on page load
app.clientside_callback(
    """
    function(n_intervals) {
        // Force dropdown z-index every few seconds to catch new elements
        const forceDropdownZIndex = () => {
            // Set high z-index on all dropdown elements
            document.querySelectorAll('#sample-selector, #gene-panel-selector').forEach(dropdown => {
                if (dropdown) {
                    dropdown.style.zIndex = '99999';
                    const parent = dropdown.closest('.glass-card');
                    if (parent) parent.style.zIndex = '1000';
                }
            });
            
            // Force all dropdown menus to highest z-index
            document.querySelectorAll('div[class*="menu"], div[class*="MenuList"], div[class*="option"], .Select-menu, .Select-menu-outer').forEach(menu => {
                menu.style.zIndex = '99999';
            });
            
            // Force React Select components
            document.querySelectorAll('.css-26l3qy-menu, .css-1pahdxg-control, .css-1hwfws3, div[class*="-menu"], div[class*="-MenuList"], div[class*="-option"]').forEach(elem => {
                elem.style.zIndex = '99999';
            });
            
            // Force all Dash dropdown components
            document.querySelectorAll('.dash-dropdown .Select-menu-outer, .dash-dropdown .Select-menu, .dash-dropdown .Select-option').forEach(elem => {
                elem.style.zIndex = '99999';
            });
        };
        
        forceDropdownZIndex();
        
        return window.dash_clientside.no_update;
    }
    """,
    Output("variant-count", "style", allow_duplicate=True),
    [Input("apply-filters-btn", "n_clicks")],
    prevent_initial_call=True
)

# Z-INDEX FIX: Enhanced clientside callback to handle dropdown z-index issues
app.clientside_callback(
    """
    function(sidebar_open) {
        setTimeout(() => {
            // SUPER AGGRESSIVE z-index forcing for ALL dropdown elements
            const forceAllDropdownsToTop = () => {
                // Target ALL possible dropdown selectors
                const allDropdownSelectors = [
                    '#sample-selector', '#gene-panel-selector',
                    '.dash-dropdown', '.Select', '.Select-control', '.Select-menu', '.Select-menu-outer',
                    '.Select-option', '.Select-input', '.Select-placeholder',
                    'div[class*="control"]', 'div[class*="menu"]', 'div[class*="MenuList"]', 
                    'div[class*="option"]', 'div[class*="placeholder"]',
                    '.css-26l3qy-menu', '.css-1pahdxg-control', '.css-1hwfws3',
                    'div[class*="-menu"]', 'div[class*="-MenuList"]', 'div[class*="-option"]',
                    'div[class*="-control"]', 'div[class*="-placeholder"]'
                ];
                
                allDropdownSelectors.forEach(selector => {
                    document.querySelectorAll(selector).forEach(elem => {
                        elem.style.zIndex = '99999 !important';
                        elem.style.position = 'relative';
                    });
                });
                
                // Ensure containers have proper z-index
                document.querySelectorAll('.glass-card').forEach(card => {
                    if (card.querySelector('#sample-selector') || card.querySelector('#gene-panel-selector')) {
                        card.style.zIndex = '1000';
                        card.style.position = 'relative';
                    }
                });
            };
            
            const sidebar = document.querySelector('.filter-sidebar');
            const overlay = document.querySelector('.sidebar-overlay');
            
            if (sidebar_open) {
                // Sidebar is open - ensure sidebar is highest but dropdowns still work
                if (sidebar) {
                    sidebar.style.zIndex = '999999';
                }
                if (overlay) {
                    overlay.style.zIndex = '999998';
                }
            }
            
            // ALWAYS force dropdowns to top regardless of sidebar state
            forceAllDropdownsToTop();
            
            // Run again multiple times to catch dynamic elements
            setTimeout(forceAllDropdownsToTop, 100);
            setTimeout(forceAllDropdownsToTop, 500);
            
        }, 50);
        
        return window.dash_clientside.no_update;
    }
    """,
    Output("sidebar-open", "data", allow_duplicate=True),
    [Input("sidebar-open", "data")],
    prevent_initial_call=True
)

# Clientside callback for filter buttons styling
app.clientside_callback(
    """
    function(active_filters) {
        setTimeout(() => {
            const filterButtons = document.querySelectorAll('[id*="preset-filter"]');
            filterButtons.forEach(button => {
                const buttonText = button.textContent.toLowerCase();
                let isActive = false;
                
                if (active_filters) {
                    if (active_filters.high_impact && buttonText.includes('high impact')) isActive = true;
                    if (active_filters.pathogenic && buttonText.includes('pathogenic')) isActive = true;
                    if (active_filters.heterozygous && buttonText.includes('heterozygous')) isActive = true;
                    if (active_filters.homozygous && buttonText.includes('homozygous')) isActive = true;
                    // NOUVEAUX : Filtres MoI
                    if (active_filters.moi_ad && buttonText.includes('autosomal dominant')) isActive = true;
                    if (active_filters.moi_ar && buttonText.includes('autosomal recessive')) isActive = true;
                    if (active_filters.moi_xl && buttonText.includes('x-linked')) isActive = true;
                }
                
                if (isActive) {
                    button.classList.remove('btn-outline-info');
                    button.classList.add('btn-info', 'active');
                    button.style.background = 'linear-gradient(45deg, #00BCD4, #0097A7)';
                    button.style.color = 'white';
                    button.style.borderColor = '#00BCD4';
                } else {
                    button.classList.remove('btn-info', 'active');
                    button.classList.add('btn-outline-info');
                    button.style.background = '';
                    button.style.color = '';
                    button.style.borderColor = '';
                }
            });
        }, 100);
        
        return window.dash_clientside.no_update;
    }
    """,
    Output("variant-count", "style"),
    [Input("active-filters", "data")]
)

# =============================================================================
# APP STARTUP AND RUN
# =============================================================================

def run_app():
    """Run the Dash application"""
    try:
        print("🚀 Starting Variant Visualizer with Gene Panel Support...")
        print("🔗 Access the app at: http://127.0.0.1:8051")
        print("📊 Database status:", "Ready" if os.path.exists(VARIANTS_PARQUET_PATH) else "No data file")
        print("🧬 Gene panels:", "Loaded" if get_available_panels() else "Empty")
        
        app.run_server(
            debug=True, 
            host='127.0.0.1', 
            port=8051,
            dev_tools_hot_reload=True,
            dev_tools_ui=True
        )
        
    except Exception as e:
        logger.error(f"Failed to start application: {e}")
        print(f"❌ Error starting app: {e}")

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 8051))
    app.run(host="0.0.0.0", port=port, debug=False)