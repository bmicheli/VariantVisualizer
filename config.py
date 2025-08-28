"""
Configuration file for Variant Visualizer
Contains all constants, paths, and configuration settings
OPTIMIZED VERSION with performance improvements
UPDATED WITH LARGER FONTS FOR BETTER READABILITY
FIXED: DROPDOWN BEHAVIOR HARMONIZATION
FIXED: ACCORDION INTERACTION ISSUES
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
# PERFORMANCE CONFIGURATION - OPTIMISATIONS PRINCIPALES
# =============================================================================

# OPTIMISATION MAJEURE: Maximum variants to display for performance
# CHANGÉ DE 1000 À 20 pour un rendu ultra-rapide
MAX_DISPLAY_VARIANTS = 20

# OPTIMISATION: Maximum variants to load at once  
# RÉDUIT de 50000 à 1000 pour améliorer la réactivité
MAX_LOAD_LIMIT = 5000

# Default chunk size for processing
DEFAULT_CHUNK_SIZE = 10000

# Cache settings - AMÉLIORÉ
CACHE_TIMEOUT = 600      # Augmenté à 10 minutes

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
    {"name": "High Impact variants", "id": "high_impact", "icon": "⚠️"},
    {"name": "Pathogenic/Likely pathogenic", "id": "pathogenic", "icon": "🔴"},
    {"name": "Heterozygous", "id": "heterozygous", "icon": "🧬"},
    {"name": "Homozygous", "id": "homozygous", "icon": "🧬"},
]

# =============================================================================
# DROPDOWN OPTIONS
# =============================================================================

CHROMOSOME_OPTIONS = [{"label": "All", "value": "all"}] + [
    {"label": f"Chr {i}", "value": str(i)} for i in range(1, 23)
] + [
    {"label": "Chr X", "value": "X"},
    {"label": "Chr Y", "value": "Y"}
]

# =============================================================================
# EXTERNAL STYLESHEETS
# =============================================================================

EXTERNAL_STYLESHEETS = [
    "https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"
]

# =============================================================================
# CSS STYLES - AVEC OPTIMISATIONS DE PERFORMANCE ET POLICES HARMONISÉES
# FIXED: ACCORDION INTERACTION ISSUES
# =============================================================================

CUSTOM_CSS = '''
body {
    background: linear-gradient(135deg, #00BCD4 0%, #4DD0E1 50%, #80E5A3 100%);
    min-height: 100vh;
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    margin: 0;
    padding: 0;
    font-size: 15px !important;  /* TAILLE DE BASE AUGMENTÉE */
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

/* DROPDOWN FIXES - SIMPLIFIED AND IDENTICAL BEHAVIOR */
.sample-selector-container,
.gene-panel-selector-container {
    position: relative !important;
    z-index: 1000 !important;
}

.sample-selector-container .Select-menu-outer,
.sample-selector-container .Select-menu,
.sample-selector-container .Select-option,
.gene-panel-selector-container .Select-menu-outer,
.gene-panel-selector-container .Select-menu,
.gene-panel-selector-container .Select-option {
    z-index: 1001 !important;
}

/* Force identical behavior for both dropdowns */
#sample-selector,
#gene-panel-selector {
    position: relative !important;
    z-index: 1000 !important;
}

#sample-selector ~ div,
#gene-panel-selector ~ div {
    position: relative !important;
    z-index: 1001 !important;
}

/* Prevent container expansion */
#sample-selector .Select-control,
#gene-panel-selector .Select-control {
    min-height: 38px !important;
    max-height: 38px !important;
    overflow: hidden !important;
}

/* Force menu positioning */
.css-26l3qy-menu,
.css-1pahdxg-control,
.css-1hwfws3,
div[class*="-menu"],
div[class*="-MenuList"],
div[class*="-option"] {
    z-index: 1002 !important;
    position: relative !important;
}

/* Main content should be lower */
.main-filters-panel {
    z-index: 100 !important;
}

/* OPTIMISATIONS DE PERFORMANCE POUR LES VARIANTS - POLICES AUGMENTÉES */
.variants-table-container {
    width: 100%;
    overflow-x: auto;
    border-radius: 12px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
    z-index: 50 !important;
    contain: layout style;
}

.variants-table {
    width: 100%;
    min-width: 1400px;
    border-collapse: collapse;
    background: white;
    table-layout: fixed;
}

.variants-table th {
    padding: 18px 15px !important;  /* Augmenté le padding */
    background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
    font-weight: 700 !important;  /* Augmenté le poids */
    font-size: 15px !important;   /* Augmenté de 13px à 15px */
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
    border-bottom: 1px solid #e9ecef;
    font-size: 14px !important;  /* Augmenté de 12px à 14px */
    line-height: 1.4 !important;  /* Améliore la lisibilité */
    /* Hauteur et padding gérés par les styles inline Python */
}

.variant-row {
    cursor: pointer;
    will-change: background-color;
    /* Hauteur gérée par les styles inline Python */
}

.variant-row:hover {
    background: rgba(0, 188, 212, 0.05) !important;
    transform: translateZ(0);
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
    font-size: 0.9em !important;  /* Augmenté de 0.85em à 0.9em */
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
    font-size: 0.85em !important;  /* Augmenté de 0.8em à 0.85em */
    padding: 3px 8px !important;   /* Augmenté le padding */
    border-radius: 8px;
    will-change: auto;
    backface-visibility: hidden;
    font-weight: 600 !important;   /* Ajouté pour plus de visibilité */
}

.gt-het { 
    background: #fff3cd; 
    color: #856404; 
    border: 1px solid #ffeaa7 !important;  /* Ajouté une bordure */
}
.gt-hom-alt { 
    background: #f8d7da; 
    color: #721c24; 
    border: 1px solid #f5c2c7 !important;  /* Ajouté une bordure */
}
.gt-hom-ref { 
    background: #d1edff; 
    color: #0c5460; 
    border: 1px solid #b6e8ff !important;  /* Ajouté une bordure */
}
.gt-missing { 
    background: #f8f9fa; 
    color: #6c757d; 
    border: 1px solid #dee2e6 !important;  /* Ajouté une bordure */
}

.detail-section {
    background: white;
    border-radius: 8px;
    padding: 18px !important;  /* Augmenté de 15px à 18px */
    margin-bottom: 15px;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    border-left: 4px solid #00BCD4;
}

.detail-section .fw-bold, .detail-section strong {
    font-size: 14px !important;  /* Augmenté pour les labels */
}

.uniform-height {
    min-height: 180px !important;
    display: flex;
    flex-direction: column;
}

.uniform-height .text-primary {
    flex-shrink: 0;
    font-size: 16px !important;  /* Augmenté pour les titres de sections */
}

.uniform-height > div:last-child {
    flex: 1;
    display: flex;
    flex-direction: column;
    justify-content: flex-start;
}

.variant-details-content {
    contain: layout style;
}

.comment-item {
    background: rgba(255, 255, 255, 0.9);
    border-radius: 8px;
    padding: 15px !important;  /* Augmenté de 12px à 15px */
    margin-bottom: 10px;
    border-left: 3px solid #00BCD4;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    font-size: 14px !important;  /* Ajouté pour les commentaires */
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

.alert-warning {
    border-left: 4px solid #f39c12;
    font-size: 0.95em !important;  /* Augmenté de 0.9em à 0.95em */
}

.alert-warning .fas {
    color: #f39c12;
}

.badge {
    will-change: auto;
    backface-visibility: hidden;
    font-size: 0.8em !important;  /* Augmenté pour les badges */
    padding: 0.4em 0.6em !important;  /* Augmenté le padding des badges */
}

/* AMÉLIORATION DES LIENS DE GÈNES */
.gene-link {
    font-size: 14px !important;  /* Augmenté pour les liens de gènes */
    font-weight: 600 !important;
}

/* AMÉLIORATION DES BOUTONS */
.btn-sm {
    font-size: 0.8em !important;  /* Augmenté légèrement */
    padding: 0.375rem 0.75rem !important;  /* Augmenté le padding */
}

/* AMÉLIORATION DE LA LISIBILITÉ GÉNÉRALE - HARMONISATION */
body, .dash-bootstrap {
    font-size: 15px !important;  /* Taille de base harmonisée */
}

/* TITRES HARMONISÉS */
h1 {
    font-size: 2.2rem !important;  /* Titre principal plus grand */
    font-weight: 700 !important;
}

h2 {
    font-size: 1.8rem !important;
    font-weight: 650 !important;
}

h3 {
    font-size: 1.5rem !important;
    font-weight: 600 !important;
}

h4 {
    font-size: 1.3rem !important;
    font-weight: 600 !important;
}

h5 {
    font-size: 1.1rem !important;
    font-weight: 600 !important;
}

h6 {
    font-size: 1.05rem !important;
    font-weight: 600 !important;
}

/* TEXTES ET LABELS */
.text-muted {
    font-size: 14px !important;
}

label, .form-label {
    font-size: 15px !important;
    font-weight: 500 !important;
}

p {
    font-size: 15px !important;
    line-height: 1.5 !important;
}

/* NAVIGATION ET MENUS */
.nav-link {
    font-size: 15px !important;
}

.navbar-brand {
    font-size: 1.4rem !important;
    font-weight: 700 !important;
}

/* BOUTONS HARMONISÉS */
.btn {
    font-size: 14px !important;
    padding: 0.6rem 1.2rem !important;
    font-weight: 500 !important;
}

.btn-sm {
    font-size: 13px !important;
    padding: 0.4rem 0.8rem !important;
}

.btn-lg {
    font-size: 16px !important;
    padding: 0.8rem 1.5rem !important;
}

/* CARDS ET CONTAINERS */
.card-title {
    font-size: 1.2rem !important;
    font-weight: 600 !important;
}

.card-text {
    font-size: 15px !important;
}

.card-body {
    padding: 1.5rem !important;
}

/* ALERTS ET NOTIFICATIONS */
.alert {
    font-size: 15px !important;
    padding: 1rem 1.25rem !important;
}

.alert h4, .alert h5 {
    font-size: 1.2rem !important;
}

/* SAMPLE SELECTOR AMÉLIORÉ */
.sample-selector-container h6 {
    font-size: 1.1rem !important;
    font-weight: 600 !important;
}

.sample-selector-container p {
    font-size: 14px !important;
}

/* QUICK FILTERS AMÉLIORÉS */
.quick-filter-btn {
    transition: all 0.2s ease;
    border-radius: 20px;
    font-size: 14px !important;  /* Augmenté */
    margin: 3px;
    padding: 0.5rem 1rem !important;  /* Augmenté */
    font-weight: 500 !important;
}

.quick-filter-btn:hover {
    transform: scale(1.05);
}

.quick-filter-btn.active {
    background: linear-gradient(45deg, #00BCD4, #0097A7) !important;
    color: white !important;
    border-color: #00BCD4 !important;
}

/* VARIANT COUNT DISPLAY */
.variant-count h5 {
    font-size: 1.4rem !important;
}

.variant-count .text-primary {
    font-size: 1.6rem !important;
}

/* AMÉLIORATION DES DROPDOWNS */
.Select-control, .css-1pahdxg-control {
    font-size: 15px !important;
    min-height: 42px !important;
    padding: 8px 12px !important;
}

.Select-option, .css-1n7v3ny-option {
    font-size: 15px !important;
    padding: 12px 16px !important;
}

.Select-placeholder, .css-1wa3eu0-placeholder {
    font-size: 15px !important;
}

.Select-value-label, .css-1uccc91-singleValue {
    font-size: 15px !important;
}

/* MODALS HARMONISÉS */
.modal-title {
    font-size: 1.4rem !important;
    font-weight: 600 !important;
}

.modal-body {
    font-size: 15px !important;
}

.modal-header {
    padding: 1.2rem 1.5rem !important;
}

.modal-body {
    padding: 1.5rem !important;
}

.modal-footer {
    padding: 1rem 1.5rem !important;
}

/* TEXTAREA ET INPUTS */
.form-control, textarea, input {
    font-size: 15px !important;
    padding: 0.6rem 0.75rem !important;
}

/* SIDEBAR AMÉLIORÉ */
.filter-sidebar h4 {
    font-size: 1.3rem !important;
}

.filter-sidebar label {
    font-size: 15px !important;
    font-weight: 600 !important;
}

.filter-sidebar .small {
    font-size: 13px !important;
}

/* RANGE SLIDER LABELS */
.rc-slider-mark-text {
    font-size: 12px !important;
}

/* SPINNER ET LOADING */
.loading-spinner h5 {
    font-size: 1.2rem !important;
}

/* DATABASE STATUS */
.db-status-container strong {
    font-size: 15px !important;
}

.db-status-container span {
    font-size: 14px !important;
}

/* TOOLTIPS */
.tooltip {
    font-size: 14px !important;
}

/* CONTAINER SPACING AMÉLIORÉ */
.container-fluid {
    padding: 25px !important;
}

/* FIXED: Styles pour empêcher les interactions indésirables avec l'accordéon */
.variant-details-content * {
    user-select: text !important;
}

.variant-details-content .detail-section {
    pointer-events: auto !important;
}

.variant-details-content a, 
.variant-details-content button, 
.variant-details-content .btn,
.variant-details-content .aa-change-toggle {
    pointer-events: auto !important;
    user-select: none !important;
}

/* FIXED: Style spécifique pour les AA change toggles */
.aa-change-toggle {
    pointer-events: auto !important;
    user-select: none !important;
    display: inline-block !important;
}

.aa-change-toggle span:last-child {
    pointer-events: auto !important;
    user-select: none !important;
}

/* FIXED: Empêcher les clics sur les détails de fermer l'accordéon parent */
.variant-details-content {
    pointer-events: auto !important;
}

.variant-details-content > * {
    pointer-events: auto !important;
}

/* FIXED: Assurer que la sélection de texte fonctionne dans les bonnes zones */
.detail-section p,
.detail-section span:not(.aa-change-toggle span),
.comment-item,
.variant-details-content .mb-2:not(.aa-change-toggle) {
    user-select: text !important;
}

/* FIXED: Éviter la propagation des événements sur les éléments interactifs */
[data-stop-propagation="true"] {
    pointer-events: auto !important;
}

/* FIXED: Style pour les collapse AA changes */
.aa-change-collapse .card-body {
    pointer-events: auto !important;
    user-select: text !important;
}

/* FIXED: Empêcher la propagation sur les éléments sélectionnables */
.variant-details-content .detail-section span,
.variant-details-content .detail-section p,
.variant-details-content .comment-item {
    pointer-events: auto !important;
}

/* FIXED: Style pour éviter les conflits d'événements */
.variant-row .variants-table td {
    position: relative;
}

.variant-details-content {
    position: relative;
    z-index: 1;
}


'''