"""
Gene Panel Management for Variant Visualizer
Handles fetching, storing, and filtering by gene panels from PanelApp UK, Australia, and internal panels
UPDATED: REMOVED GREEN GENE FILTERING AND FIXED CONFIDENCE LEVEL PARSING
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
                            'gene_confidence': 'GREEN',  # Internal panels are all green
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
        """Fetch ALL panels from PanelApp UK with pagination support"""
        logger.info("Fetching ALL panels from PanelApp UK with pagination...")
        panels_data = []
        
        try:
            # Start with first page to get total count
            url = f"{PANELAPP_UK_BASE_URL}/panels/"
            page = 1
            total_panels_fetched = 0
            
            while url:
                logger.info(f"Fetching UK panels page {page}...")
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                panels_list = data.get('results', [])
                total_count = data.get('count', 0)
                next_url = data.get('next')
                
                if page == 1:
                    logger.info(f"Total panels available in PanelApp UK: {total_count}")
                
                logger.info(f"Processing {len(panels_list)} panels from page {page}")
                
                for i, panel in enumerate(panels_list):
                    panel_id = str(panel['id'])
                    panel_name = panel['name']
                    panel_version = panel.get('version', '1.0')
                    
                    total_panels_fetched += 1
                    logger.info(f"Fetching UK panel {total_panels_fetched}/{total_count}: {panel_name}")
                    
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
                        
                        # Debug: Log first few genes to understand structure
                        if total_panels_fetched == 1 and genes:
                            logger.info(f"Sample UK gene structure: {genes[0]}")
                            # Log ALL fields to understand the structure
                            logger.info(f"Sample UK gene ALL FIELDS: {json.dumps(genes[0], indent=2)}")
                        
                        # Format panel name: "Panel Name Version X (Y genes)"
                        formatted_name = f"{panel_name} Version {panel_version} ({gene_count} genes)"
                        
                        for gene in genes:
                            # DEBUG: Log the first few genes completely to understand structure
                            if total_panels_fetched == 1 and len(panels_data) < 5:
                                logger.info(f"DEBUG UK Gene: {json.dumps(gene, indent=2)}")
                            
                            # Parse gene data - try different possible field structures
                            gene_data = gene.get('gene_data', {})
                            gene_symbol = gene_data.get('gene_symbol', '') or gene_data.get('hgnc_symbol', '') or gene.get('gene_symbol', '')
                            
                            # Parse confidence level - ENHANCED DEBUG VERSION
                            confidence_raw = None
                            confidence_field = None
                            
                            # Try all possible confidence fields
                            possible_fields = ['confidence_level', 'confidence', 'level_of_confidence', 'rating', 'gene_confidence']
                            for field in possible_fields:
                                if field in gene:
                                    confidence_raw = gene[field]
                                    confidence_field = field
                                    break
                                elif field in gene_data:
                                    confidence_raw = gene_data[field]
                                    confidence_field = f"gene_data.{field}"
                                    break
                            
                            # Log what we found for debugging
                            if total_panels_fetched == 1 and len(panels_data) < 5:
                                logger.info(f"Gene {gene_symbol}: confidence_field='{confidence_field}', confidence_raw='{confidence_raw}'")
                                logger.info(f"Available gene fields: {list(gene.keys())}")
                                logger.info(f"Available gene_data fields: {list(gene_data.keys())}")
                            
                            # Normalize confidence values
                            if confidence_raw:
                                confidence_upper = str(confidence_raw).upper()
                                if 'GREEN' in confidence_upper or confidence_upper == '3':
                                    confidence_level = 'GREEN'
                                elif 'AMBER' in confidence_upper or confidence_upper == '2':
                                    confidence_level = 'AMBER'
                                elif 'RED' in confidence_upper or confidence_upper == '1':
                                    confidence_level = 'RED'
                                else:
                                    confidence_level = 'UNKNOWN'
                                    if total_panels_fetched <= 3:  # Only log for first few panels to avoid spam
                                        logger.warning(f"Unknown UK confidence level: '{confidence_raw}' (type: {type(confidence_raw)}) for gene {gene_symbol}")
                            else:
                                confidence_level = 'UNKNOWN'
                                if total_panels_fetched <= 3:
                                    logger.warning(f"No confidence field found for UK gene {gene_symbol}")
                            
                            if gene_symbol:  # Only add if we have a gene symbol
                                panels_data.append({
                                    'panel_id': f"uk_{panel_id}",
                                    'panel_name': formatted_name,
                                    'source': 'panelapp_uk',
                                    'gene_symbol': gene_symbol,
                                    'gene_confidence': confidence_level,
                                    'panel_version': panel_version,
                                    'last_updated': datetime.now().isoformat(),
                                    'panel_url': f"https://panelapp.genomicsengland.co.uk/panels/{panel_id}/"
                                })
                        
                        # Log confidence level distribution for first few panels
                        if total_panels_fetched <= 3:
                            recent_entries = [p for p in panels_data if p['panel_id'] == f"uk_{panel_id}"]
                            conf_counts = {}
                            for entry in recent_entries:
                                conf = entry['gene_confidence']
                                conf_counts[conf] = conf_counts.get(conf, 0) + 1
                            logger.info(f"UK Panel {panel_name} confidence distribution: {conf_counts}")
                        
                        # Add delay to be respectful to API
                        time.sleep(0.3)  # Reduced delay for faster processing
                        
                    except Exception as e:
                        logger.error(f"Error fetching UK panel {panel_id}: {e}")
                        continue
                
                # Move to next page
                url = next_url
                page += 1
                
                # Log progress every 10 pages
                if page % 10 == 0:
                    logger.info(f"Processed {total_panels_fetched} panels so far...")
                    
        except Exception as e:
            logger.error(f"Error fetching PanelApp UK panels: {e}")
        
        logger.info(f"Fetched {len(panels_data)} gene entries from {total_panels_fetched} PanelApp UK panels")
        return panels_data
    
    def fetch_panelapp_au_panels(self, max_panels=None):
        """Fetch ALL panels from PanelApp Australia with pagination support"""
        logger.info("Fetching ALL panels from PanelApp Australia with pagination...")
        panels_data = []
        
        try:
            # Start with first page to get total count
            url = f"{PANELAPP_AU_BASE_URL}/panels/"
            page = 1
            total_panels_fetched = 0
            
            while url:
                logger.info(f"Fetching AU panels page {page}...")
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                panels_list = data.get('results', [])
                total_count = data.get('count', 0)
                next_url = data.get('next')
                
                if page == 1:
                    logger.info(f"Total panels available in PanelApp Australia: {total_count}")
                
                logger.info(f"Processing {len(panels_list)} panels from page {page}")
                
                for i, panel in enumerate(panels_list):
                    panel_id = str(panel['id'])
                    panel_name = panel['name']
                    panel_version = panel.get('version', '1.0')
                    
                    total_panels_fetched += 1
                    logger.info(f"Fetching AU panel {total_panels_fetched}/{total_count}: {panel_name}")
                    
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
                        
                        # Debug: Log first few genes to understand structure
                        if total_panels_fetched == 1 and genes:
                            logger.info(f"Sample AU gene structure: {genes[0]}")
                            # Log ALL fields to understand the structure
                            logger.info(f"Sample AU gene ALL FIELDS: {json.dumps(genes[0], indent=2)}")
                        
                        # Format panel name: "Panel Name Version X (Y genes)"
                        formatted_name = f"{panel_name} Version {panel_version} ({gene_count} genes)"
                        
                        for gene in genes:
                            # DEBUG: Log the first few genes completely to understand structure
                            if total_panels_fetched == 1 and len(panels_data) < 5:
                                logger.info(f"DEBUG AU Gene: {json.dumps(gene, indent=2)}")
                            
                            # Parse gene data - try different possible field structures
                            gene_data = gene.get('gene_data', {})
                            gene_symbol = gene_data.get('gene_symbol', '') or gene_data.get('hgnc_symbol', '') or gene.get('gene_symbol', '')
                            
                            # Parse confidence level - ENHANCED DEBUG VERSION
                            confidence_raw = None
                            confidence_field = None
                            
                            # Try all possible confidence fields
                            possible_fields = ['confidence_level', 'confidence', 'level_of_confidence', 'rating', 'gene_confidence']
                            for field in possible_fields:
                                if field in gene:
                                    confidence_raw = gene[field]
                                    confidence_field = field
                                    break
                                elif field in gene_data:
                                    confidence_raw = gene_data[field]
                                    confidence_field = f"gene_data.{field}"
                                    break
                            
                            # Log what we found for debugging
                            if total_panels_fetched == 1 and len(panels_data) < 5:
                                logger.info(f"Gene {gene_symbol}: confidence_field='{confidence_field}', confidence_raw='{confidence_raw}'")
                                logger.info(f"Available gene fields: {list(gene.keys())}")
                                logger.info(f"Available gene_data fields: {list(gene_data.keys())}")
                            
                            # Normalize confidence values
                            if confidence_raw:
                                confidence_upper = str(confidence_raw).upper()
                                if 'GREEN' in confidence_upper or confidence_upper == '3':
                                    confidence_level = 'GREEN'
                                elif 'AMBER' in confidence_upper or confidence_upper == '2':
                                    confidence_level = 'AMBER'
                                elif 'RED' in confidence_upper or confidence_upper == '1':
                                    confidence_level = 'RED'
                                else:
                                    confidence_level = 'UNKNOWN'
                                    if total_panels_fetched <= 3:  # Only log for first few panels to avoid spam
                                        logger.warning(f"Unknown AU confidence level: '{confidence_raw}' (type: {type(confidence_raw)}) for gene {gene_symbol}")
                            else:
                                confidence_level = 'UNKNOWN'
                                if total_panels_fetched <= 3:
                                    logger.warning(f"No confidence field found for AU gene {gene_symbol}")
                            
                            if gene_symbol:  # Only add if we have a gene symbol
                                panels_data.append({
                                    'panel_id': f"au_{panel_id}",
                                    'panel_name': formatted_name,
                                    'source': 'panelapp_au',
                                    'gene_symbol': gene_symbol,
                                    'gene_confidence': confidence_level,
                                    'panel_version': panel_version,
                                    'last_updated': datetime.now().isoformat(),
                                    'panel_url': f"https://panelapp-aus.org/panels/{panel_id}/"
                                })
                        
                        # Log confidence level distribution for first few panels
                        if total_panels_fetched <= 3:
                            recent_entries = [p for p in panels_data if p['panel_id'] == f"au_{panel_id}"]
                            conf_counts = {}
                            for entry in recent_entries:
                                conf = entry['gene_confidence']
                                conf_counts[conf] = conf_counts.get(conf, 0) + 1
                            logger.info(f"AU Panel {panel_name} confidence distribution: {conf_counts}")
                        
                        # Add delay to be respectful to API
                        time.sleep(0.3)  # Reduced delay for faster processing
                        
                    except Exception as e:
                        logger.error(f"Error fetching AU panel {panel_id}: {e}")
                        continue
                
                # Move to next page
                url = next_url
                page += 1
                
                # Log progress every 10 pages  
                if page % 10 == 0:
                    logger.info(f"Processed {total_panels_fetched} panels so far...")
                    
        except Exception as e:
            logger.error(f"Error fetching PanelApp Australia panels: {e}")
        
        logger.info(f"Fetched {len(panels_data)} gene entries from {total_panels_fetched} PanelApp Australia panels")
        return panels_data
    
    def update_all_panels(self, max_panels_per_source=None):
        """Update all panels from external sources and reload internal panels"""
        logger.info("Starting full panel update...")
        
        all_panels_data = []
        
        # Always reload internal panels first
        logger.info("Reloading internal panels...")
        internal_panels_data = []
        self._load_internal_panels_to_list(internal_panels_data)
        all_panels_data.extend(internal_panels_data)
        
        # Fetch from PanelApp UK - GET ALL PANELS
        uk_panels = self.fetch_panelapp_uk_panels(max_panels_per_source)
        all_panels_data.extend(uk_panels)
        
        # Fetch from PanelApp Australia - GET ALL PANELS  
        au_panels = self.fetch_panelapp_au_panels(max_panels_per_source)
        all_panels_data.extend(au_panels)
        
        # Create new DataFrame
        if all_panels_data:
            self.panels_df = pl.DataFrame(all_panels_data)
            self.save_panels()
            
            # Log overall confidence distribution - FIXED VERSION
            if len(self.panels_df) > 0:
                # Use manual counting instead of value_counts to avoid index issues
                confidence_counts = {}
                for _, row in self.panels_df.to_pandas().iterrows():
                    conf = row['gene_confidence']
                    confidence_counts[conf] = confidence_counts.get(conf, 0) + 1
                logger.info(f"Overall confidence distribution: {confidence_counts}")
            
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
                            'gene_confidence': 'GREEN',  # Internal panels are all green
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
    
    def get_genes_for_panels(self, panel_ids):
        """Get gene list for selected panels - SIMPLIFIED (no green gene filtering)"""
        if not panel_ids or len(self.panels_df) == 0:
            return []
        
        try:
            genes = (
                self.panels_df
                .filter(pl.col('panel_id').is_in(panel_ids))
                .select('gene_symbol')
                .unique()
                .to_pandas()['gene_symbol']
                .tolist()
            )
            
            # Remove empty strings and duplicates
            genes = list(set([g.strip() for g in genes if g and g.strip()]))
            
            logger.info(f"Found {len(genes)} unique genes for {len(panel_ids)} panels")
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
                
                # Get gene count and confidence distribution
                panel_genes = self.panels_df.filter(pl.col('panel_id') == panel_id)
                total_genes = len(panel_genes.select('gene_symbol').unique())
                
                # Get confidence distribution - FIXED VERSION
                confidence_counts = {}
                for _, row in panel_genes.to_pandas().iterrows():
                    conf = row['gene_confidence']
                    confidence_counts[conf] = confidence_counts.get(conf, 0) + 1
                
                green_genes = confidence_counts.get('GREEN', 0)
                amber_genes = confidence_counts.get('AMBER', 0) 
                red_genes = confidence_counts.get('RED', 0)
                unknown_genes = confidence_counts.get('UNKNOWN', 0)
                
                panel_info['gene_count'] = total_genes
                panel_info['green_gene_count'] = green_genes
                panel_info['amber_gene_count'] = amber_genes
                panel_info['red_gene_count'] = red_genes
                panel_info['unknown_gene_count'] = unknown_genes
                panel_info['confidence_distribution'] = confidence_counts  # Use the fixed dict
                
                # Debug log for this specific panel
                logger.info(f"Panel {panel_id} info: Total={total_genes}, Green={green_genes}, Amber={amber_genes}, Red={red_genes}, Unknown={unknown_genes}")
                
                return panel_info
            
        except Exception as e:
            logger.error(f"Error getting panel info: {e}")
        
        return None
    
    def get_panel_count(self):
        """Get total number of unique panels"""
        if len(self.panels_df) == 0:
            return 0
        
        return len(self.panels_df.select('panel_id').unique())
    
    def should_update(self, days_threshold=2):
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
            panel_manager.update_all_panels()  # REMOVED max_panels limit
        except Exception as e:
            logger.error(f"Error updating panels: {e}")

def get_available_panels():
    """Get available panels for UI"""
    return panel_manager.get_available_panels()

def get_genes_for_panels_optimized(panel_ids):
    """OPTIMIZED: Get genes for selected panels - MUCH FASTER for large panels"""
    if not panel_ids or len(panel_manager.panels_df) == 0:
        return []
    
    try:
        # OPTIMIZATION: Use Polars native operations for speed
        genes = (
            panel_manager.panels_df
            .filter(pl.col('panel_id').is_in(panel_ids))
            .select('gene_symbol')
            .unique()
            .to_series()
            .to_list()
        )
        
        # OPTIMIZATION: Use set operations for deduplication (faster than list comprehension)
        genes_set = {g.strip() for g in genes if g and g.strip()}
        unique_genes = list(genes_set)
        
        logger.info(f"OPTIMIZED: Found {len(unique_genes)} unique genes for {len(panel_ids)} panels")
        return unique_genes
        
    except Exception as e:
        logger.error(f"Error getting genes for panels: {e}")
        return []

def get_genes_for_panels(panel_ids):
    """Get genes for selected panels - uses optimized version"""
    return get_genes_for_panels_optimized(panel_ids)

def get_panel_info(panel_id):
    """Get panel information"""
    return panel_manager.get_panel_info(panel_id)

def force_update_panels():
    """Force update all panels - GET ALL PANELS"""
    panel_manager.update_all_panels()  # REMOVED max_panels limit