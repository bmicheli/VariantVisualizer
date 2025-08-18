"""
Main application file for Variant Visualizer
Contains the Dash app initialization, layout, and all callbacks
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

# =============================================================================
# APP INITIALIZATION
# =============================================================================

# Initialize Parquet database
print("🔧 Initializing Parquet database...")
init_parquet_database()
print("✅ Parquet database structure initialized.")

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=EXTERNAL_STYLESHEETS, suppress_callback_exceptions=True)

# Set app title
app.title = "Variant Visualizer - ClinVar Focus"

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
        
        # Database Status Panel
        dbc.Card([
            dbc.CardBody([
                html.H6([html.I(className="fas fa-database me-2"), "Database Status"], className="mb-2 text-primary"),
                html.Div(id="database-status-info")
            ], style={"padding": "15px"})
        ], className="db-status-container mb-3"),
        
        # Sample Selector Panel
        create_sample_selector(),
        
        # Main Filters Panel
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
    create_stats_modal(),
    
    # Data stores
    dcc.Store(id="filtered-variants"),
    dcc.Store(id="selected-variant-id"),
    dcc.Store(id="active-filters", data={}),
    dcc.Store(id="sidebar-open", data=False),
    dcc.Store(id="available-samples", data=[]),
    
    # Download component
    dcc.Download(id="download-csv")
    
], style={"minHeight": "100vh", "background": "linear-gradient(135deg, #00BCD4 0%, #4DD0E1 50%, #80E5A3 100%)"})

# =============================================================================
# CALLBACKS
# =============================================================================

# Database status callback
@app.callback(
    Output("database-status-info", "children"),
    [Input("refresh-data-btn", "n_clicks")]
)
def update_database_status(refresh_clicks):
    """Update database status information"""
    if refresh_clicks:
        db._invalidate_cache()
    return create_database_status_display()

# Sample options callback
@app.callback(
    [Output("sample-selector", "options"), Output("available-samples", "data")],
    [Input("refresh-data-btn", "n_clicks")]
)
def update_sample_options(refresh_clicks):
    """Update sample selector options when data is refreshed"""
    if refresh_clicks:
        db._invalidate_cache()
    
    samples = get_available_samples()
    options = [{"label": sample, "value": sample} for sample in samples]
    return options, samples

# Show all samples callback
@app.callback(
    Output("sample-selector", "value"),
    [Input("show-all-samples", "n_clicks"), Input("clear-samples", "n_clicks")],
    [State("available-samples", "data")],
    prevent_initial_call=True
)
def handle_sample_selection(show_all_clicks, clear_clicks, available_samples):
    """Handle sample selection buttons"""
    ctx = callback_context
    if not ctx.triggered:
        return dash.no_update
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == "show-all-samples" and show_all_clicks and available_samples:
        return available_samples
    elif button_id == "clear-samples" and clear_clicks:
        return []
    
    return dash.no_update

# Sidebar toggle callback
@app.callback(
    [Output("filter-sidebar", "className"), Output("sidebar-overlay", "className"), Output("sidebar-open", "data")],
    [Input("more-filters-btn", "n_clicks"), Input("close-sidebar-btn", "n_clicks"), Input("sidebar-overlay", "n_clicks")],
    [State("sidebar-open", "data")],
    prevent_initial_call=True
)
def toggle_sidebar(toggle_clicks, close_clicks, overlay_clicks, is_open):
    """Toggle sidebar visibility"""
    ctx = callback_context
    if not ctx.triggered:
        return "filter-sidebar", "sidebar-overlay", False
    
    new_state = not is_open
    sidebar_class = "filter-sidebar open" if new_state else "filter-sidebar"
    overlay_class = "sidebar-overlay open" if new_state else "sidebar-overlay"
    return sidebar_class, overlay_class, new_state

# Main variants display callback
@app.callback(
    [Output("variants-display", "children"), Output("variant-count", "children"), Output("filtered-variants", "data")],
    [Input("search-input", "value"), Input("genotype-filter", "value"), Input("chromosome-filter", "value"),
     Input("active-filters", "data"), Input("sample-selector", "value"), Input("refresh-data-btn", "n_clicks"),
     Input("quality-filter", "value"), Input("vaf-range-slider", "value")]
)
def update_variants_display(search_term, genotype_filter, chromosome_filter, active_filters, 
                          selected_samples, refresh_clicks, quality_filter, vaf_range):
    """Update the variants display with optimized filtering"""
    
    if refresh_clicks:
        db._invalidate_cache()
    
    # Early return if no samples selected
    if not selected_samples:
        no_selection = create_no_selection_display()
        count_display = create_variant_count_display(0, 0, 0)
        return no_selection, count_display, []
    
    try:
        # Load variants with basic sample filtering at the database level
        df = db.load_variants_lazy(samples=selected_samples, limit=MAX_LOAD_LIMIT)
        
        if is_dataframe_empty(df):
            empty_display = create_beautiful_variant_display(df)
            empty_count = create_variant_count_display(0, 0, len(selected_samples))
            return empty_display, empty_count, []
        
        original_count = len(df)
        
        # Apply additional filters efficiently
        df = apply_filters(
            df, 
            search_term=search_term,
            genotype_filter=genotype_filter,
            chromosome_filter=chromosome_filter,
            active_filters=active_filters,
            selected_samples=selected_samples
        )
        
        # Apply advanced filters
        if isinstance(df, pd.DataFrame):
            df = pl.from_pandas(df)
        
        if quality_filter is not None and quality_filter > 0:
            df = df.filter(pl.col('QUAL') >= quality_filter)
        
        if vaf_range and len(vaf_range) == 2:
            df = df.filter((pl.col('VAF') >= vaf_range[0]) & (pl.col('VAF') <= vaf_range[1]))
        
        display = create_beautiful_variant_display(df)
        count_display = create_variant_count_display(len(df), original_count, len(selected_samples))
        
        # Convert to records for storage (limit to prevent memory issues)
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

# Clear filters callback
@app.callback(
    [Output("genotype-filter", "value"), Output("chromosome-filter", "value"), 
     Output("search-input", "value"), Output("active-filters", "data", allow_duplicate=True),
     Output("quality-filter", "value"), Output("vaf-range-slider", "value")],
    [Input("clear-filters", "n_clicks")],
    prevent_initial_call=True
)
def clear_all_filters(n_clicks):
    """Clear all filters"""
    if n_clicks:
        return "all", "all", "", {}, None, [0, 1]
    return dash.no_update

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
            # Keep modal open on error
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

# Statistics modal callback
@app.callback(
    [Output("stats-modal", "is_open"), Output("stats-content", "children")],
    [Input("stats-btn", "n_clicks"), Input("close-stats-modal", "n_clicks")],
    [State("stats-modal", "is_open"), State("filtered-variants", "data")]
)
def handle_stats_modal(stats_clicks, close_clicks, is_open, filtered_data):
    """Handle statistics modal"""
    ctx = callback_context
    if not ctx.triggered:
        return False, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == "stats-btn" and stats_clicks:
        # Create statistics from filtered data
        if filtered_data:
            df = pd.DataFrame(filtered_data)
            stats = create_summary_stats(df)
        else:
            stats = db.get_database_stats()
        
        stats_display = create_stats_display(stats)
        return True, stats_display
    
    if trigger_id == "close-stats-modal":
        return False, dash.no_update
    
    return is_open, dash.no_update

# Export callback
@app.callback(
    Output("download-csv", "data"),
    [Input("export-btn", "n_clicks")],
    [State("filtered-variants", "data")],
    prevent_initial_call=True
)
def export_variants(n_clicks, filtered_data):
    """Export filtered variants to CSV"""
    if n_clicks and filtered_data:
        df = pd.DataFrame(filtered_data)
        export_columns = [
            'CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT', 'VAF', 'gene', 
            'consequence', 'clinvar_sig', 'clinvar_id', 'clinvar_disease',
            'gnomad_af', 'cadd_score', 'review_status'
        ]
        available_columns = [col for col in export_columns if col in df.columns]
        export_df = df[available_columns].copy()
        
        filename = f"variants_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        
        return dcc.send_data_frame(
            export_df.to_csv, 
            filename, 
            index=False
        )
    return dash.no_update

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
                    if (active_filters.rare && buttonText.includes('rare')) isActive = true;
                    if (active_filters.high_impact && buttonText.includes('high impact')) isActive = true;
                    if (active_filters.reviewed && buttonText.includes('reviewed')) isActive = true;
                    if (active_filters.pending && buttonText.includes('pending')) isActive = true;
                    if (active_filters.clinvar_annotated && buttonText.includes('clinvar annotated')) isActive = true;
                    if (active_filters.pathogenic && buttonText.includes('pathogenic')) isActive = true;
                    if (active_filters.benign && buttonText.includes('benign')) isActive = true;
                    if (active_filters.vus && buttonText.includes('vus')) isActive = true;
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

# Error handling callback
@app.callback(
    Output("variants-display", "children", allow_duplicate=True),
    [Input("refresh-data-btn", "n_clicks")],
    prevent_initial_call=True
)
def handle_refresh_errors(n_clicks):
    """Handle errors during data refresh"""
    if n_clicks:
        try:
            # Test database connection
            stats = db.get_database_stats()
            if stats['total_variants'] == 0:
                return create_error_component(
                    "No variants found in database. Please check your data files."
                )
        except Exception as e:
            logger.error(f"Database error during refresh: {e}")
            return create_error_component(
                f"Database connection error: {str(e)}"
            )
    
    return dash.no_update

# Real-time search callback (with debouncing)
@app.callback(
    Output("search-input", "valid"),
    [Input("search-input", "value")],
    prevent_initial_call=True
)
def validate_search_input(search_value):
    """Validate search input in real-time"""
    if search_value and len(search_value) < 2:
        return False
    return True

# Performance monitoring callback
@app.callback(
    Output("database-status-info", "children", allow_duplicate=True),
    [Input("variants-display", "children")],
    prevent_initial_call=True
)
def monitor_performance(variants_display):
    """Monitor app performance and update status"""
    try:
        # Get current stats
        stats = db.get_database_stats()
        
        # Add performance info
        if hasattr(db, '_last_query_time'):
            performance_info = html.Div([
                create_database_status_display(),
                html.Small(
                    f"Last query: {db._last_query_time:.2f}s", 
                    className="text-muted mt-1"
                )
            ])
            return performance_info
        
        return create_database_status_display()
        
    except Exception as e:
        logger.error(f"Performance monitoring error: {e}")
        return create_database_status_display()

# =============================================================================
# APP STARTUP AND RUN
# =============================================================================

def run_app():
    """Run the Dash application"""
    try:
        print("🚀 Starting Variant Visualizer...")
        print("🔗 Access the app at: http://127.0.0.1:8050")
        print("📊 Database status:", "Ready" if os.path.exists(VARIANTS_PARQUET_PATH) else "No data file")
        
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

