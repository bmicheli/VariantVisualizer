"""
Configuration file for Variant Visualizer
Contains all constants, paths, and configuration settings
"""

import os
import logging

# =============================================================================
# FILE PATHS AND DIRECTORIES
# =============================================================================

DATA_DIR = 'data'
VARIANTS_PARQUET_PATH = os.path.join(DATA_DIR, 'variants.parquet')
COMMENTS_PARQUET_PATH = os.path.join(DATA_DIR, 'comments.parquet')
SAMPLE_INDEX_PATH = os.path.join(DATA_DIR, 'sample_index.parquet')

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# =============================================================================
# UI CONFIGURATION
# =============================================================================

# Color mappings for consequences
CONSEQUENCE_COLORS = {
    'frameshift_variant': 'danger',
    'stop_gained': 'danger',
    'stopgain': 'danger', 
    'stop_lost': 'danger',
    'start_lost': 'danger',
    'splice_acceptor_variant': 'danger',
    'splice_donor_variant': 'danger',
    'missense_variant': 'warning',
    'nonsynonymous_SNV': 'warning',
    'synonymous_variant': 'success',
    'synonymous_SNV': 'success',
    'synonymous': 'success',
    'frameshift_deletion': 'danger',
    'frameshift_insertion': 'danger',
    'variant': 'secondary'
}

# Color mappings for ClinVar classifications
CLINVAR_COLORS = {
    'Pathogenic': 'danger',
    'Likely pathogenic': 'warning', 
    'VUS': 'secondary',
    'Uncertain significance': 'secondary',
    'Likely benign': 'info',
    'Benign': 'success',
    'Conflicting': 'dark',
    'Other': 'light',
    'Not provided': 'light',
    'Drug response': 'primary',
    'Association': 'info',
    'Protective': 'success',
    'Risk factor': 'warning'
}

# Color mappings for review status
STATUS_COLORS = {
    'Reviewed': 'success',
    'Pending': 'warning'
}

# =============================================================================
# FILTER PRESETS
# =============================================================================

PRESET_FILTERS = [
    {"name": "Rare variants (AF < 0.1%)", "id": "rare", "icon": "💎"},
    {"name": "High Impact variants", "id": "high_impact", "icon": "⚠️"},
    {"name": "ClinVar annotated", "id": "clinvar_annotated", "icon": "📋"},
    {"name": "Pathogenic/Likely pathogenic", "id": "pathogenic", "icon": "🔴"},
    {"name": "VUS", "id": "vus", "icon": "❓"},
    {"name": "Benign/Likely benign", "id": "benign", "icon": "🟢"},
    {"name": "Reviewed variants", "id": "reviewed", "icon": "✅"},
    {"name": "Pending review", "id": "pending", "icon": "⏳"},
]

# =============================================================================
# DROPDOWN OPTIONS
# =============================================================================

GENOTYPE_OPTIONS = [
    {"label": "All", "value": "all"},
    {"label": "Heterozygous (0/1)", "value": "het"},
    {"label": "Homozygous Alt (1/1)", "value": "hom_alt"},
    {"label": "Homozygous Ref (0/0)", "value": "hom_ref"}
]

CHROMOSOME_OPTIONS = [{"label": "All", "value": "all"}] + [
    {"label": f"Chr {i}", "value": str(i)} for i in range(1, 23)
] + [
    {"label": "Chr X", "value": "X"},
    {"label": "Chr Y", "value": "Y"}
]

# =============================================================================
# DATABASE SETTINGS
# =============================================================================

# Default chunk size for processing
DEFAULT_CHUNK_SIZE = 10000

# Maximum variants to display for performance
MAX_DISPLAY_VARIANTS = 1000

# Maximum variants to load at once
MAX_LOAD_LIMIT = 50000

# Cache settings
CACHE_TIMEOUT = 300  # 5 minutes

# =============================================================================
# EXTERNAL STYLESHEETS
# =============================================================================

EXTERNAL_STYLESHEETS = [
    "https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"
]

# =============================================================================
# CSS STYLES
# =============================================================================

CUSTOM_CSS = '''
body {
    background: linear-gradient(135deg, #00BCD4 0%, #4DD0E1 50%, #80E5A3 100%);
    min-height: 100vh;
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    margin: 0;
    padding: 0;
}

.glass-card {
    background: rgba(255, 255, 255, 0.95) !important;
    backdrop-filter: blur(10px);
    border-radius: 15px !important;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1) !important;
    border: 1px solid rgba(255, 255, 255, 0.2) !important;
}

/* SIDEBAR Z-INDEX FIXES */
.filter-sidebar {
    position: fixed !important;
    top: 0 !important;
    left: -400px !important;
    width: 380px !important;
    height: 100vh !important;
    background: rgba(255, 255, 255, 0.98) !important;
    backdrop-filter: blur(15px) !important;
    box-shadow: 2px 0 20px rgba(0, 0, 0, 0.15) !important;
    transition: left 0.3s ease !important;
    z-index: 99999 !important;
    overflow-y: auto !important;
    border-radius: 0 15px 15px 0 !important;
}

.filter-sidebar.open {
    left: 0 !important;
    z-index: 99999 !important;
}

.sidebar-overlay {
    position: fixed !important;
    top: 0 !important;
    left: 0 !important;
    width: 100vw !important;
    height: 100vh !important;
    background: rgba(0, 0, 0, 0.5) !important;
    opacity: 0 !important;
    visibility: hidden !important;
    transition: all 0.3s ease !important;
    z-index: 99998 !important;
}

.sidebar-overlay.open {
    opacity: 1 !important;
    visibility: visible !important;
    z-index: 99998 !important;
}

/* DROPDOWN Z-INDEX FIXES */
.sample-selector-container {
    position: relative !important;
    z-index: 1000 !important;
}

.sample-selector-container .Select-menu-outer,
.sample-selector-container .Select-menu,
.sample-selector-container .Select-option {
    z-index: 1001 !important;
}

/* Dash specific dropdown fixes */
.dash-dropdown {
    z-index: 1001 !important;
}

.dash-dropdown .Select-menu-outer {
    z-index: 1002 !important;
}

.dash-dropdown .Select-menu {
    z-index: 1003 !important;
}

.dash-dropdown .Select-option {
    z-index: 1004 !important;
}

/* React-Select fixes */
.css-26l3qy-menu,
.css-1pahdxg-control,
.css-1hwfws3 {
    z-index: 1005 !important;
}

div[class*="-menu"] {
    z-index: 1005 !important;
}

div[class*="-MenuList"] {
    z-index: 1006 !important;
}

div[class*="-option"] {
    z-index: 1007 !important;
}

/* Main content should be lower */
.main-filters-panel {
    z-index: 100 !important;
}

.variants-table-container {
    width: 100%;
    overflow-x: auto;
    border-radius: 12px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
    z-index: 50 !important;
}

.variants-table {
    width: 100%;
    min-width: 1400px;
    border-collapse: collapse;
    background: white;
}

.variants-table th {
    padding: 15px 12px;
    background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
    font-weight: 600;
    font-size: 13px;
    text-align: left;
    border-bottom: 2px solid #dee2e6;
    white-space: nowrap;
    position: sticky;
    top: 0;
    z-index: 10;
}

.sortable-header {
    transition: background-color 0.2s ease;
    cursor: pointer;
    user-select: none;
}

.sortable-header:hover {
    background: linear-gradient(135deg, #e9ecef 0%, #dee2e6 100%) !important;
}

.variants-table td {
    padding: 12px;
    border-bottom: 1px solid #e9ecef;
    vertical-align: middle;
    font-size: 12px;
}

.variant-row {
    cursor: pointer;
}

.variant-row:hover {
    background: rgba(0, 188, 212, 0.05) !important;
}

.main-filters-panel {
    background: rgba(255, 255, 255, 0.95) !important;
    backdrop-filter: blur(10px);
    border-radius: 15px !important;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1) !important;
    border: 1px solid rgba(255, 255, 255, 0.2) !important;
}

.quick-filter-btn {
    transition: all 0.2s ease;
    border-radius: 20px;
    font-size: 0.85em;
    margin: 2px;
}

.quick-filter-btn:hover {
    transform: scale(1.05);
}

.quick-filter-btn.active {
    background: linear-gradient(45deg, #00BCD4, #0097A7) !important;
    color: white !important;
    border-color: #00BCD4 !important;
}

.sample-selector-container, .db-status-container {
    background: rgba(255, 255, 255, 0.95) !important;
    backdrop-filter: blur(10px);
    border-radius: 12px !important;
    border: 1px solid rgba(255, 255, 255, 0.2) !important;
    box-shadow: 0 2px 15px rgba(0, 0, 0, 0.08) !important;
}

.genotype-badge {
    font-family: monospace;
    font-size: 0.8em;
    padding: 2px 6px;
    border-radius: 8px;
}

.gt-het { background: #fff3cd; color: #856404; }
.gt-hom-alt { background: #f8d7da; color: #721c24; }
.gt-hom-ref { background: #d1edff; color: #0c5460; }
.gt-missing { background: #f8f9fa; color: #6c757d; }

.detail-section {
    background: white;
    border-radius: 8px;
    padding: 15px;
    margin-bottom: 15px;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    border-left: 4px solid #00BCD4;
}

/* Uniform height fix for the six boxes */
.uniform-height {
    min-height: 180px !important;
    display: flex;
    flex-direction: column;
}

.uniform-height .text-primary {
    flex-shrink: 0;
}

.uniform-height > div:last-child {
    flex: 1;
    display: flex;
    flex-direction: column;
    justify-content: flex-start;
}

.comment-item {
    background: rgba(255, 255, 255, 0.9);
    border-radius: 8px;
    padding: 12px;
    margin-bottom: 10px;
    border-left: 3px solid #00BCD4;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.loading-spinner {
    display: flex;
    justify-content: center;
    align-items: center;
    height: 200px;
}

.error-message {
    color: #dc3545;
    background: rgba(220, 53, 69, 0.1);
    border: 1px solid rgba(220, 53, 69, 0.2);
    border-radius: 8px;
    padding: 15px;
    margin: 10px 0;
}

.success-message {
    color: #198754;
    background: rgba(25, 135, 84, 0.1);
    border: 1px solid rgba(25, 135, 84, 0.2);
    border-radius: 8px;
    padding: 15px;
    margin: 10px 0;
}
'''
