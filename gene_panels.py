"""
Gene Panel Management for Variant Visualizer
Handles fetching, storing, and filtering by gene panels from PanelApp UK, Australia, and internal panels
UPDATED WITH CORRECT AU API AND SIMPLIFIED DISPLAY
"""

import requests
import pandas as pd
import polars as pl
import json
import os
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
from config import DATA_DIR, logger
import time
from pathlib import Path

# Configuration
PANELAPP_UK_BASE_URL = "https://panelapp.genomicsengland.co.uk/api/v1"
PANELAPP_AU_BASE_URL = "https://panelapp-aus.org/api/v1"
GENE_PANELS_DB_PATH = os.path.join(DATA_DIR, 'gene_panels.parquet')
PANEL_METADATA_PATH = os.path.join(DATA_DIR, 'panel_metadata.json')
INTERNAL_PANELS_DIR = os.path.join(DATA_DIR, 'internal_panels')

class GenePanelManager:
    """Manages gene panels from multiple sources"""
    
    def __init__(self):
        self.panels_df = None
        self.metadata = {}
        self.load_cached_panels()
    
    def load_cached_panels(self):
        """Load cached gene panels from disk"""
        try:
            if os.path.exists(GENE_PANELS_DB_PATH):
                self.panels_df = pl.read_parquet(GENE_PANELS_DB_PATH)
                logger.info(f"Loaded {len(self.panels_df)} gene panel entries from cache")
                
                if os.path.exists(PANEL_METADATA_PATH):
                    with open(PANEL_METADATA_PATH, 'r') as f:
                        self.metadata = json.load(f)
            else:
                logger.info("No cached gene panels found, loading internal panels...")
                self.panels_df = self._create_empty_panels_df()
                self._load_internal_panels()
                
        except Exception as e:
            logger.error(f"Error loading cached panels: {e}")
            self.panels_df = self._create_empty_panels_df()
            self._load_internal_panels()
    
    def _create_empty_panels_df(self):
        """Create empty DataFrame with proper schema"""
        schema = {
            'panel_id': pl.Utf8,
            'panel_name': pl.Utf8,
            'source': pl.Utf8,      # 'panelapp_uk', 'panelapp_au', 'internal'
            'gene_symbol': pl.Utf8,
            'gene_confidence': pl.Utf8,  # 'GREEN', 'AMBER', 'RED'
            'panel_version': pl.Utf8,
            'last_updated': pl.Utf8,
            'panel_url': pl.Utf8
        }
        return pl.DataFrame(schema=schema)
    
    def _load_internal_panels(self):
        """Load internal panels from text files in internal_panels directory"""
        if not os.path.exists(INTERNAL_PANELS_DIR):
            logger.warning(f"Internal panels directory not found: {INTERNAL_PANELS_DIR}")
            return
        
        logger.info("Loading internal panels from text files...")
        panels_data = []
        
        try:
            # Get all .txt files in the internal_panels directory
            panel_files = list(Path(INTERNAL_PANELS_DIR).glob("*.txt"))
            logger.info(f"Found {len(panel_files)} internal panel files")
            
            for panel_file in panel_files:
                panel_name = self._parse_panel_name(panel_file.stem)
                panel_id = f"internal_{panel_file.stem.lower()}"
                
                logger.info(f"Loading panel: {panel_name}")
                
                try:
                    with open(panel_file, 'r', encoding='utf-8') as f:
                        genes = []
                        for line in f:
                            gene = line.strip()
                            if gene and not gene.startswith('#'):  # Skip empty lines and comments
                                genes.append(gene)
                    
                    logger.info(f"  Found {len(genes)} genes in {panel_name}")
                    
                    # Create panel entries for each gene
                    for gene_symbol in genes:
                        panels_data.append({
                            'panel_id': panel_id,
                            'panel_name': panel_name,
                            'source': 'internal',
                            'gene_symbol': gene_symbol.strip(),
                            'gene_confidence': 'GREEN',
                            'panel_version': self._extract_version(panel_file.stem),
                            'last_updated': datetime.now().isoformat(),
                            'panel_url': ''
                        })
                        
                except Exception as e:
                    logger.error(f"Error reading panel file {panel_file}: {e}")
                    continue
            
            if panels_data:
                self.panels_df = pl.DataFrame(panels_data)
                self.save_panels()
                logger.info(f"Successfully loaded {len(panels_data)} gene entries from {len(panel_files)} internal panels")
            else:
                logger.warning("No internal panel data loaded")
                
        except Exception as e:
            logger.error(f"Error loading internal panels: {e}")
    
    def _parse_panel_name(self, filename):
        """Parse panel name from filename with simplified format"""
        # Remove version info and clean up name
        name = filename.replace('_', ' ')
        
        # Extract version and number
        import re
        
        # Extract version (like _v1, _v2, etc.)
        version_match = re.search(r'(\d+)_v(\d+)$', filename)
        if version_match:
            gene_count = version_match.group(1)
            version = version_match.group(2)
            # Remove the count and version from the name
            name = re.sub(r'_\d+_v\d+$', '', filename).replace('_', ' ')
            
            # Capitalize words
            words = name.split()
            capitalized_words = [word.capitalize() for word in words]
            clean_name = ' '.join(capitalized_words)
            
            # Format: "Panel Name Version X (Y genes)"
            return f"{clean_name} Version {version} ({gene_count} genes)"
        else:
            # Fallback for files without version pattern
            words = name.split()
            capitalized_words = [word.capitalize() if not word.isdigit() else word for word in words]
            return ' '.join(capitalized_words)
    
    def _extract_version(self, filename):
        """Extract version from filename"""
        import re
        version_match = re.search(r'_v(\d+)$', filename)
        if version_match:
            return f"v{version_match.group(1)}"
        return "v1"
    
    def fetch_panelapp_uk_panels(self, max_panels=None):
        """Fetch panels from PanelApp UK"""
        logger.info("Fetching panels from PanelApp UK...")
        panels_data = []
        
        try:
            # Get list of panels
            response = requests.get(f"{PANELAPP_UK_BASE_URL}/panels/", timeout=30)
            response.raise_for_status()
            panels_list = response.json()['results']
            
            if max_panels:
                panels_list = panels_list[:max_panels]
            
            for i, panel in enumerate(panels_list):
                panel_id = str(panel['id'])
                panel_name = panel['name']
                panel_version = panel.get('version', '1.0')
                
                logger.info(f"Fetching UK panel {i+1}/{len(panels_list)}: {panel_name}")
                
                try:
                    # Get detailed panel info with genes
                    detail_response = requests.get(
                        f"{PANELAPP_UK_BASE_URL}/panels/{panel_id}/",
                        timeout=30
                    )
                    detail_response.raise_for_status()
                    panel_detail = detail_response.json()
                    
                    # Extract genes and count them
                    genes = panel_detail.get('genes', [])
                    gene_count = len(genes)
                    
                    # Format panel name: "Panel Name Version X (Y genes)"
                    formatted_name = f"{panel_name} Version {panel_version} ({gene_count} genes)"
                    
                    for gene in genes:
                        gene_data = gene.get('gene_data', {})
                        panels_data.append({
                            'panel_id': f"uk_{panel_id}",
                            'panel_name': formatted_name,
                            'source': 'panelapp_uk',
                            'gene_symbol': gene_data.get('gene_symbol', ''),
                            'gene_confidence': gene.get('confidence_level', 'UNKNOWN'),
                            'panel_version': panel_version,
                            'last_updated': datetime.now().isoformat(),
                            'panel_url': f"https://panelapp.genomicsengland.co.uk/panels/{panel_id}/"
                        })
                    
                    # Add delay to be respectful to API
                    time.sleep(0.5)
                    
                except Exception as e:
                    logger.error(f"Error fetching UK panel {panel_id}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Error fetching PanelApp UK panels: {e}")
        
        logger.info(f"Fetched {len(panels_data)} gene entries from PanelApp UK")
        return panels_data
    
    def fetch_panelapp_au_panels(self, max_panels=None):
        """Fetch panels from PanelApp Australia"""
        logger.info("Fetching panels from PanelApp Australia...")
        panels_data = []
        
        try:
            # Get list of panels
            response = requests.get(f"{PANELAPP_AU_BASE_URL}/panels/", timeout=30)
            response.raise_for_status()
            panels_list = response.json()['results']
            
            if max_panels:
                panels_list = panels_list[:max_panels]
            
            for i, panel in enumerate(panels_list):
                panel_id = str(panel['id'])
                panel_name = panel['name']
                panel_version = panel.get('version', '1.0')
                
                logger.info(f"Fetching AU panel {i+1}/{len(panels_list)}: {panel_name}")
                
                try:
                    # Get detailed panel info with genes
                    detail_response = requests.get(
                        f"{PANELAPP_AU_BASE_URL}/panels/{panel_id}/",
                        timeout=30
                    )
                    detail_response.raise_for_status()
                    panel_detail = detail_response.json()
                    
                    # Extract genes and count them
                    genes = panel_detail.get('genes', [])
                    gene_count = len(genes)
                    
                    # Format panel name: "Panel Name Version X (Y genes)"
                    formatted_name = f"{panel_name} Version {panel_version} ({gene_count} genes)"
                    
                    for gene in genes:
                        gene_data = gene.get('gene_data', {})
                        panels_data.append({
                            'panel_id': f"au_{panel_id}",
                            'panel_name': formatted_name,
                            'source': 'panelapp_au',
                            'gene_symbol': gene_data.get('gene_symbol', ''),
                            'gene_confidence': gene.get('confidence_level', 'UNKNOWN'),
                            'panel_version': panel_version,
                            'last_updated': datetime.now().isoformat(),
                            'panel_url': f"https://panelapp-aus.org/panels/{panel_id}/"
                        })
                    
                    # Add delay to be respectful to API
                    time.sleep(0.5)
                    
                except Exception as e:
                    logger.error(f"Error fetching AU panel {panel_id}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Error fetching PanelApp Australia panels: {e}")
        
        logger.info(f"Fetched {len(panels_data)} gene entries from PanelApp Australia")
        return panels_data
    
    def update_all_panels(self, max_panels_per_source=50):
        """Update all panels from external sources and reload internal panels"""
        logger.info("Starting full panel update...")
        
        all_panels_data = []
        
        # Always reload internal panels first
        logger.info("Reloading internal panels...")
        internal_panels_data = []
        self._load_internal_panels_to_list(internal_panels_data)
        all_panels_data.extend(internal_panels_data)
        
        # Fetch from PanelApp UK
        uk_panels = self.fetch_panelapp_uk_panels(max_panels_per_source)
        all_panels_data.extend(uk_panels)
        
        # Fetch from PanelApp Australia
        au_panels = self.fetch_panelapp_au_panels(max_panels_per_source)
        all_panels_data.extend(au_panels)
        
        # Create new DataFrame
        if all_panels_data:
            self.panels_df = pl.DataFrame(all_panels_data)
            self.save_panels()
            
            # Update metadata
            self.metadata = {
                'last_update': datetime.now().isoformat(),
                'total_panels': self.get_panel_count(),
                'sources': {
                    'panelapp_uk': len(uk_panels),
                    'panelapp_au': len(au_panels),
                    'internal': len(internal_panels_data)
                }
            }
            
            with open(PANEL_METADATA_PATH, 'w') as f:
                json.dump(self.metadata, f, indent=2)
            
            logger.info(f"Panel update complete. Total gene entries: {len(all_panels_data)}")
        else:
            logger.warning("No panels fetched during update")
    
    def _load_internal_panels_to_list(self, panels_data_list):
        """Load internal panels into the provided list"""
        if not os.path.exists(INTERNAL_PANELS_DIR):
            return
        
        try:
            panel_files = list(Path(INTERNAL_PANELS_DIR).glob("*.txt"))
            
            for panel_file in panel_files:
                panel_name = self._parse_panel_name(panel_file.stem)
                panel_id = f"internal_{panel_file.stem.lower()}"
                
                try:
                    with open(panel_file, 'r', encoding='utf-8') as f:
                        genes = []
                        for line in f:
                            gene = line.strip()
                            if gene and not gene.startswith('#'):
                                genes.append(gene)
                    
                    for gene_symbol in genes:
                        panels_data_list.append({
                            'panel_id': panel_id,
                            'panel_name': panel_name,
                            'source': 'internal',
                            'gene_symbol': gene_symbol.strip(),
                            'gene_confidence': 'GREEN',
                            'panel_version': self._extract_version(panel_file.stem),
                            'last_updated': datetime.now().isoformat(),
                            'panel_url': ''
                        })
                        
                except Exception as e:
                    logger.error(f"Error reading panel file {panel_file}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Error loading internal panels to list: {e}")
    
    def save_panels(self):
        """Save panels to disk"""
        try:
            if len(self.panels_df) > 0:
                self.panels_df.write_parquet(GENE_PANELS_DB_PATH)
                logger.info(f"Saved {len(self.panels_df)} gene panel entries to disk")
        except Exception as e:
            logger.error(f"Error saving panels: {e}")
    
    def get_available_panels(self):
        """Get list of available panels for UI with simplified format"""
        if len(self.panels_df) == 0:
            return []
        
        try:
            panels = (
                self.panels_df
                .select(['panel_id', 'panel_name', 'source'])
                .unique()
                .sort(['source', 'panel_name'])
                .to_pandas()
            )
            
            panel_options = []
            for _, panel in panels.iterrows():
                source_icon = {
                    'panelapp_uk': '🇬🇧',
                    'panelapp_au': '🇦🇺', 
                    'internal': '🏠'
                }.get(panel['source'], '📋')
                
                # Simplified format: Icon + Panel Name
                label = f"{source_icon} {panel['panel_name']}"
                panel_options.append({
                    'label': label,
                    'value': panel['panel_id']
                })
            
            return panel_options
            
        except Exception as e:
            logger.error(f"Error getting available panels: {e}")
            return []
    
    def get_genes_for_panels(self, panel_ids, green_genes_only=False):
        """Get gene list for selected panels with optional green gene filtering"""
        if not panel_ids or len(self.panels_df) == 0:
            return []
        
        try:
            df_filter = pl.col('panel_id').is_in(panel_ids)
            
            # Add green gene filter if requested
            if green_genes_only:
                df_filter = df_filter & (pl.col('gene_confidence') == 'GREEN')
            
            genes = (
                self.panels_df
                .filter(df_filter)
                .select('gene_symbol')
                .unique()
                .to_pandas()['gene_symbol']
                .tolist()
            )
            
            # Remove empty strings and duplicates
            genes = list(set([g.strip() for g in genes if g and g.strip()]))
            
            filter_text = " (green genes only)" if green_genes_only else ""
            logger.info(f"Found {len(genes)} unique genes for {len(panel_ids)} panels{filter_text}")
            return genes
            
        except Exception as e:
            logger.error(f"Error getting genes for panels: {e}")
            return []
    
    def get_panel_info(self, panel_id):
        """Get detailed information about a panel"""
        if len(self.panels_df) == 0:
            return None
        
        try:
            panel_data = (
                self.panels_df
                .filter(pl.col('panel_id') == panel_id)
                .select(['panel_name', 'source', 'panel_url', 'panel_version'])
                .unique()
                .to_pandas()
            )
            
            if len(panel_data) > 0:
                panel_info = panel_data.iloc[0].to_dict()
                
                # Get gene count and green gene count
                panel_genes = self.panels_df.filter(pl.col('panel_id') == panel_id)
                total_genes = len(panel_genes.select('gene_symbol').unique())
                green_genes = len(panel_genes.filter(pl.col('gene_confidence') == 'GREEN').select('gene_symbol').unique())
                
                panel_info['gene_count'] = total_genes
                panel_info['green_gene_count'] = green_genes
                
                return panel_info
            
        except Exception as e:
            logger.error(f"Error getting panel info: {e}")
        
        return None
    
    def get_panel_count(self):
        """Get total number of unique panels"""
        if len(self.panels_df) == 0:
            return 0
        
        return len(self.panels_df.select('panel_id').unique())
    
    def should_update(self, days_threshold=7):
        """Check if panels should be updated based on last update time"""
        if not self.metadata.get('last_update'):
            return True
        
        try:
            last_update = datetime.fromisoformat(self.metadata['last_update'])
            days_since_update = (datetime.now() - last_update).days
            return days_since_update >= days_threshold
        except:
            return True
    
    def search_panels(self, search_term):
        """Search panels by name"""
        if not search_term or len(self.panels_df) == 0:
            return []
        
        try:
            search_lower = search_term.lower()
            matching_panels = (
                self.panels_df
                .filter(pl.col('panel_name').str.to_lowercase().str.contains(search_lower))
                .select(['panel_id', 'panel_name', 'source'])
                .unique()
                .to_pandas()
            )
            
            return matching_panels.to_dict('records')
            
        except Exception as e:
            logger.error(f"Error searching panels: {e}")
            return []

# Global panel manager instance
panel_manager = GenePanelManager()

def init_gene_panels():
    """Initialize gene panel system"""
    os.makedirs(DATA_DIR, exist_ok=True)
    
    # Load internal panels if no cached panels exist
    if panel_manager.get_panel_count() == 0:
        logger.info("Loading internal panels for the first time...")
        panel_manager._load_internal_panels()

def update_panels_if_needed():
    """Update panels if they're older than threshold"""
    if panel_manager.should_update():
        logger.info("Panels need updating...")
        try:
            panel_manager.update_all_panels(max_panels_per_source=30)
        except Exception as e:
            logger.error(f"Error updating panels: {e}")

def get_available_panels():
    """Get available panels for UI"""
    return panel_manager.get_available_panels()

def get_genes_for_panels(panel_ids, green_genes_only=False):
    """Get genes for selected panels with optional green gene filtering"""
    return panel_manager.get_genes_for_panels(panel_ids, green_genes_only)

def get_panel_info(panel_id):
    """Get panel information"""
    return panel_manager.get_panel_info(panel_id)

def force_update_panels():
    """Force update all panels"""
    panel_manager.update_all_panels()