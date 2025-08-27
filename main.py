"""
Main application file for Variant Visualizer
Contains the Dash app initialization, layout, and all callbacks
OPTIMIZED VERSION with performance improvements
UPDATED WITH GENE PANEL FILTERING - REMOVED GREEN GENE FILTERING
UPDATED WITH NEW GNOMAD AND CGEN FREQUENCY DISPLAY
"""

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
from gene_panels import (
    init_gene_panels, update_panels_if_needed, get_available_panels, 
    get_genes_for_panels, get_panel_info, force_update_panels
)

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
    dcc.Store(id="panel-genes", data=[])
    
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

# SIMPLIFIED: Gene panel selection callback (removed green gene filtering)
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
     Input("apply-filters-btn", "n_clicks")],  # Add apply button as trigger
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

# UPDATED: Main variants display callback - OPTIMIZED for large gene panels
@app.callback(
    [Output("variants-display", "children"), Output("variant-count", "children"), Output("filtered-variants", "data")],
    [Input("apply-filters-btn", "n_clicks"),  # Primary trigger: apply button
     Input("sample-selector", "value"),       # Also update when samples change
     Input("active-filters", "data"),         # Also update when preset filters change
     Input("panel-genes", "data")],           # Also update when gene panels change
    [State("search-input", "value"), 
     State("vaf-range-slider", "value"),
     State("selected-gene-panels", "data")],     # Get current filter values
    prevent_initial_call=False
)
def update_variants_display_optimized(apply_clicks, selected_samples, active_filters, panel_genes, 
                                    search_term, vaf_range, selected_panels):
    """SUPER OPTIMIZED variants display for large gene panels"""
    
    # Early return if no samples selected
    if not selected_samples:
        no_selection = create_no_selection_display()
        count_display = create_variant_count_display(0, 0, 0)
        return no_selection, count_display, []
    
    try:
        # PERFORMANCE: Log start time
        import time
        start_time = time.time()
        
        # OPTIMISATION 1: Use gene panel filter at database level if possible
        if panel_genes and len(panel_genes) > 100:  # Only for large panels
            logger.info(f"PERFORMANCE: Large panel detected ({len(panel_genes)} genes), using optimized filtering")
            
            # Load variants normally - db is the variants database, not panels
            df = db.load_variants_lazy(samples=selected_samples, limit=MAX_LOAD_LIMIT)
            
            if is_dataframe_empty(df):
                empty_display = create_beautiful_variant_display(df)
                empty_count = create_variant_count_display(0, 0, len(selected_samples))
                return empty_display, empty_count, []
            
            original_count = len(df)
            logger.info(f"PERFORMANCE: Loaded {original_count} variants in {time.time() - start_time:.2f}s")
            
            # OPTIMISATION 2: Pre-convert gene panel to set for O(1) lookup - SUPER OPTIMIZED
            panel_genes_set = set(panel_genes)
            logger.info(f"PERFORMANCE: Initial gene set created ({len(panel_genes_set)} entries) in {time.time() - start_time:.2f}s")
            
            # OPTIMISATION 3: Pre-compute gene mapping ONCE and create reverse lookup
            gene_mapping = load_gene_mapping()
            if gene_mapping:
                mapping_start = time.time()
                
                # Create reverse mapping: gene_name -> gene_id (much faster)
                reverse_mapping = {mapped_name.upper(): gene_id for gene_id, mapped_name in gene_mapping.items()}
                
                # Add gene IDs to the set using vectorized lookup
                for gene_name in panel_genes:
                    gene_name_upper = gene_name.upper()
                    if gene_name_upper in reverse_mapping:
                        panel_genes_set.add(reverse_mapping[gene_name_upper])
                
                logger.info(f"PERFORMANCE: Gene mapping applied in {time.time() - mapping_start:.2f}s")
            
            logger.info(f"PERFORMANCE: Final gene set prepared ({len(panel_genes_set)} entries) in {time.time() - start_time:.2f}s")
            
            # OPTIMISATION 4: Apply gene filter with vectorized operations - SUPER FAST
            if isinstance(df, pd.DataFrame):
                df = pl.from_pandas(df)
            
            # Convert set to list once for Polars filter
            filter_start = time.time()
            gene_list = list(panel_genes_set)
            
            # Use highly optimized Polars filter
            df = df.filter(pl.col('gene').is_in(gene_list))
            logger.info(f"PERFORMANCE: Gene panel filter applied ({len(df)} variants) in {time.time() - filter_start:.2f}s")
            
        else:
            # Standard loading for small panels or no panels
            df = db.load_variants_lazy(samples=selected_samples, limit=MAX_LOAD_LIMIT)
            
            if is_dataframe_empty(df):
                empty_display = create_beautiful_variant_display(df)
                empty_count = create_variant_count_display(0, 0, len(selected_samples))
                return empty_display, empty_count, []
            
            original_count = len(df)
            logger.info(f"Standard loading: {original_count} variants")
            
            # Apply gene panel filtering for smaller panels using the original method
            if panel_genes and selected_panels:
                logger.info(f"Applying standard gene panel filter with {len(panel_genes)} genes")
                
                if isinstance(df, pd.DataFrame):
                    df = pl.from_pandas(df)
                
                # Create conditions for gene filtering
                gene_conditions = []
                
                # Direct gene symbol match
                gene_conditions.append(pl.col('gene').is_in(panel_genes))
                
                # Also try to match gene names after ID conversion (if gene mapping exists)
                gene_mapping = load_gene_mapping()
                if gene_mapping:
                    # Get gene IDs that correspond to our gene names
                    matching_ids = []
                    for gene_name in panel_genes:
                        for gene_id, mapped_name in gene_mapping.items():
                            if mapped_name.upper() == gene_name.upper():
                                matching_ids.append(gene_id)
                    
                    if matching_ids:
                        gene_conditions.append(pl.col('gene').is_in(matching_ids))
                
                # Apply gene filtering
                if gene_conditions:
                    combined_gene_filter = gene_conditions[0]
                    for condition in gene_conditions[1:]:
                        combined_gene_filter = combined_gene_filter | condition
                    
                    df = df.filter(combined_gene_filter)
                    logger.info(f"After gene panel filter: {len(df)} variants")
        
        # OPTIMISATION 4: Apply remaining filters efficiently
        filter_start = time.time()
        
        # Apply basic filters (optimized)
        df = apply_filters_optimized(
            df, 
            search_term=search_term,
            genotype_filter=None,
            chromosome_filter=None,
            active_filters=active_filters,
            selected_samples=selected_samples
        )
        
        logger.info(f"PERFORMANCE: Basic filters applied in {time.time() - filter_start:.2f}s ({len(df)} variants)")
        
        # Apply VAF filter
        if isinstance(df, pd.DataFrame):
            df = pl.from_pandas(df)
        
        if vaf_range and len(vaf_range) == 2:
            before_vaf = len(df)
            df = df.filter((pl.col('VAF') >= vaf_range[0]) & (pl.col('VAF') <= vaf_range[1]))
            logger.info(f"VAF filter ({vaf_range}): {before_vaf} -> {len(df)}")
        
        total_time = time.time() - start_time
        logger.info(f"PERFORMANCE: Total filtering time: {total_time:.2f}s for {len(df)} final variants")
        
        display = create_beautiful_variant_display(df)
        count_display = create_variant_count_display(len(df), original_count, len(selected_samples))
        
        # OPTIMISATION 5: Store only first MAX_DISPLAY_VARIANTS for details
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

# Reset all filters callback - MAIN PANEL RESET BUTTON
@app.callback(
    [Output("search-input", "value", allow_duplicate=True), 
     Output("active-filters", "data", allow_duplicate=True),
     Output("vaf-range-slider", "value", allow_duplicate=True)],
    [Input("reset-all-btn", "n_clicks")],
    prevent_initial_call=True
)
def reset_all_filters(n_clicks):
    """Reset all filters from the main panel reset button"""
    if n_clicks:
        return "", {}, [0, 1]
    return dash.no_update

# Updated Clear filters callback to include apply functionality
@app.callback(
    [Output("search-input", "value", allow_duplicate=True), 
     Output("active-filters", "data", allow_duplicate=True),
     Output("vaf-range-slider", "value", allow_duplicate=True),
     Output("filter-sidebar", "className", allow_duplicate=True), 
     Output("sidebar-overlay", "className", allow_duplicate=True), 
     Output("sidebar-open", "data", allow_duplicate=True)],
    [Input("clear-filters", "n_clicks")],
    prevent_initial_call=True
)
def clear_all_filters_sidebar(n_clicks):
    """Clear filters from sidebar clear button and close sidebar"""
    if n_clicks:
        return "", {}, [0, 1], "filter-sidebar", "sidebar-overlay", False
    return dash.no_update

# Variant row expansion callback
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
        print("🔗 Access the app at: http://127.0.0.1:8050")
        print("📊 Database status:", "Ready" if os.path.exists(VARIANTS_PARQUET_PATH) else "No data file")
        print("🧬 Gene panels:", "Loaded" if get_available_panels() else "Empty")
        
        app.run_server(
            debug=True, 
            host='127.0.0.1', 
            port=8050,
            dev_tools_hot_reload=True,
            dev_tools_ui=True
        )
        
    except Exception as e:
        logger.error(f"Failed to start application: {e}")
        print(f"❌ Error starting app: {e}")

if __name__ == "__main__":
    app.run(debug=True)