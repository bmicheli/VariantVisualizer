"""
Database operations for Variant Visualizer
Handles Parquet file operations, caching, and data management
OPTIMIZED VERSION with performance improvements
UPDATED WITH NEW GNOMAD AND CGEN FREQUENCY FIELDS
"""

import polars as pl
import pandas as pd
import os
from datetime import datetime
from functools import lru_cache
from typing import Dict, List, Optional, Any
from config import *

class OptimizedParquetDB:
    """Optimized Parquet database manager with caching and lazy loading"""
    
    def __init__(self):
        self.variants_df = None
        self.sample_index = None
        self.comments_df = None
        self._last_modified = None
        self._cache = {}
        
    @lru_cache(maxsize=128)
    def get_available_samples(self):
        """Get cached list of available samples"""
        try:
            if os.path.exists(SAMPLE_INDEX_PATH):
                sample_df = pl.read_parquet(SAMPLE_INDEX_PATH)
                return sorted(sample_df['SAMPLE'].to_list())
            elif os.path.exists(VARIANTS_PARQUET_PATH):
                df = pl.read_parquet(VARIANTS_PARQUET_PATH, columns=['SAMPLE'])
                return sorted(df['SAMPLE'].unique().to_list())
        except Exception as e:
            logger.error(f"Error loading samples: {e}")
        return []
    
    def _check_file_modified(self):
        """Check if parquet files have been modified"""
        if not os.path.exists(VARIANTS_PARQUET_PATH):
            return False
            
        current_modified = os.path.getmtime(VARIANTS_PARQUET_PATH)
        if self._last_modified != current_modified:
            self._last_modified = current_modified
            self._invalidate_cache()
            return True
        return False
    
    def _invalidate_cache(self):
        """Invalidate internal caches"""
        self.variants_df = None
        self.sample_index = None
        self._cache.clear()
        self.get_available_samples.cache_clear()
    
    def load_variants_lazy(self, samples=None, chromosomes=None, limit=None):
        """Optimized lazy load variants with aggressive filtering"""
        try:
            if not os.path.exists(VARIANTS_PARQUET_PATH):
                return pl.DataFrame()
            
            # OPTIMISATION: Use more aggressive effective limit
            effective_limit = min(limit or MAX_LOAD_LIMIT, MAX_LOAD_LIMIT)
            
            # Build filter conditions
            filters = []
            if samples:
                filters.append(pl.col('SAMPLE').is_in(samples))
            if chromosomes:
                chrom_list = [str(c).replace('chr', '') for c in chromosomes]
                filters.append(pl.col('CHROM').is_in(chrom_list))
            
            # Read with lazy evaluation and early limiting
            lazy_df = pl.scan_parquet(VARIANTS_PARQUET_PATH)
            
            # Apply most selective filters first
            if filters:
                for f in filters:
                    lazy_df = lazy_df.filter(f)
            
            # OPTIMISATION: Apply limit very early
            lazy_df = lazy_df.limit(effective_limit)
            
            # Collect results
            df = lazy_df.collect()
            
            # Ensure max_gnomad_af is calculated if missing (backward compatibility)
            if len(df) > 0:
                df = self._ensure_max_gnomad_af(df)
            
            # OPTIMISATION: Conditional comment loading only if necessary
            if len(df) <= MAX_DISPLAY_VARIANTS:
                if os.path.exists(COMMENTS_PARQUET_PATH):
                    try:
                        comments_df = pl.read_parquet(COMMENTS_PARQUET_PATH)
                        if len(comments_df) > 0:
                            comment_counts = (
                                comments_df
                                .group_by(['variant_key', 'sample_id'])
                                .agg(pl.len().alias('comment_count'))
                            )
                            df = df.join(
                                comment_counts, 
                                left_on=['variant_key', 'SAMPLE'], 
                                right_on=['variant_key', 'sample_id'],
                                how='left'
                            ).with_columns(
                                pl.col('comment_count').fill_null(0)
                            )
                    except Exception as e:
                        logger.warning(f"Could not load comment counts: {e}")
                        df = df.with_columns(pl.lit(0).alias('comment_count'))
                else:
                    df = df.with_columns(pl.lit(0).alias('comment_count'))
            else:
                # If too many variants, don't load comments for now
                df = df.with_columns(pl.lit(0).alias('comment_count'))
            
            return df
            
        except Exception as e:
            logger.error(f"Error loading variants: {e}")
            return pl.DataFrame()
    
    def _ensure_max_gnomad_af(self, df):
        """Ensure max_gnomad_af column exists and is calculated correctly"""
        try:
            # Check if max_gnomad_af column exists
            if 'max_gnomad_af' not in df.columns:
                logger.info("Calculating max_gnomad_af from population-specific frequencies")
                
                # Get available gnomAD population columns (excluding ASJ and FIN)
                gnomad_pop_columns = []
                for col in ['gnomad_af_afr', 'gnomad_af_amr', 'gnomad_af_eas', 'gnomad_af_nfe', 'gnomad_af_sas']:
                    if col in df.columns:
                        gnomad_pop_columns.append(col)
                
                if gnomad_pop_columns:
                    # Calculate max across available population columns
                    df = df.with_columns([
                        pl.max_horizontal([pl.col(col) for col in gnomad_pop_columns]).alias('max_gnomad_af')
                    ])
                else:
                    # Fallback to general gnomad_af or 0
                    if 'gnomad_af' in df.columns:
                        df = df.with_columns([
                            pl.col('gnomad_af').alias('max_gnomad_af')
                        ])
                    else:
                        df = df.with_columns([
                            pl.lit(0.0).alias('max_gnomad_af')
                        ])
            
            return df
            
        except Exception as e:
            logger.error(f"Error calculating max_gnomad_af: {e}")
            # Add zero column as fallback
            return df.with_columns([pl.lit(0.0).alias('max_gnomad_af')])
    
    def get_database_stats(self):
        """Get database statistics efficiently"""
        cache_key = "db_stats"
        if cache_key in self._cache and not self._check_file_modified():
            return self._cache[cache_key]
        
        try:
            if not os.path.exists(VARIANTS_PARQUET_PATH):
                stats = {
                    'total_variants': 0,
                    'total_samples': 0,
                    'reviewed_variants': 0,
                    'pending_variants': 0,
                    'clinvar_annotated': 0,
                    'chromosomes': []
                }
            else:
                # Use lazy scanning for efficiency
                lazy_df = pl.scan_parquet(VARIANTS_PARQUET_PATH)
                
                # Get basic counts
                total_variants = lazy_df.select(pl.len()).collect().item()
                total_samples = lazy_df.select(pl.col('SAMPLE').n_unique()).collect().item()
                
                # Get review status counts
                status_counts = (
                    lazy_df
                    .group_by('review_status')
                    .agg(pl.len().alias('count'))
                    .collect()
                )
                
                reviewed = status_counts.filter(pl.col('review_status') == 'Reviewed')['count'].sum()
                pending = status_counts.filter(pl.col('review_status') == 'Pending')['count'].sum()
                
                # Get ClinVar annotated count
                clinvar_annotated = (
                    lazy_df
                    .filter(pl.col('clinvar_sig').is_not_null())
                    .select(pl.len())
                    .collect().item()
                )
                
                # Get chromosomes
                chromosomes = lazy_df.select(pl.col('CHROM').unique()).collect()['CHROM'].to_list()
                
                stats = {
                    'total_variants': total_variants,
                    'total_samples': total_samples,
                    'reviewed_variants': reviewed or 0,
                    'pending_variants': pending or total_variants,
                    'clinvar_annotated': clinvar_annotated or 0,
                    'chromosomes': sorted(chromosomes)
                }
            
            self._cache[cache_key] = stats
            return stats
            
        except Exception as e:
            logger.error(f"Error getting database stats: {e}")
            return {
                'total_variants': 0,
                'total_samples': 0,
                'reviewed_variants': 0,
                'pending_variants': 0,
                'clinvar_annotated': 0,
                'chromosomes': []
            }

    def search_variants(self, search_term: str, limit: int = 1000):
        """Search variants by gene, consequence, or sample"""
        try:
            if not os.path.exists(VARIANTS_PARQUET_PATH):
                return pl.DataFrame()
            
            search_lower = search_term.lower()
            lazy_df = pl.scan_parquet(VARIANTS_PARQUET_PATH)
            
            search_conditions = [
                pl.col('gene').str.to_lowercase().str.contains(search_lower),
                pl.col('consequence').str.to_lowercase().str.contains(search_lower),
                pl.col('SAMPLE').str.to_lowercase().str.contains(search_lower),
                pl.col('CHROM').str.to_lowercase().str.contains(search_lower)
            ]
            
            # Combine search conditions with OR
            combined_search = search_conditions[0]
            for condition in search_conditions[1:]:
                combined_search = combined_search | condition
            
            df = lazy_df.filter(combined_search).limit(limit).collect()
            
            # Ensure max_gnomad_af is available
            if len(df) > 0:
                df = self._ensure_max_gnomad_af(df)
            
            return df
            
        except Exception as e:
            logger.error(f"Error searching variants: {e}")
            return pl.DataFrame()

# Global database instance
db = OptimizedParquetDB()

def init_parquet_database():
    """Initialize Parquet database structure"""
    os.makedirs(DATA_DIR, exist_ok=True)
    
    if not os.path.exists(COMMENTS_PARQUET_PATH):
        comments_schema = {
            'id': pl.Series([], dtype=pl.Int64),
            'variant_key': pl.Series([], dtype=pl.Utf8),
            'sample_id': pl.Series([], dtype=pl.Utf8),
            'user_name': pl.Series([], dtype=pl.Utf8),
            'comment_text': pl.Series([], dtype=pl.Utf8),
            'timestamp': pl.Series([], dtype=pl.Utf8)
        }
        empty_comments_df = pl.DataFrame(comments_schema)
        empty_comments_df.write_parquet(COMMENTS_PARQUET_PATH)

def load_variants_from_parquet():
    """Load variants with optimized filtering and caching"""
    return db.load_variants_lazy()

def get_available_samples():
    """Get list of unique sample IDs from Parquet database"""
    return db.get_available_samples()

def add_variant_comment(variant_key, sample_id, user_name, comment_text):
    """Add comment to variant in Parquet database"""
    try:
        if os.path.exists(COMMENTS_PARQUET_PATH):
            comments_df = pl.read_parquet(COMMENTS_PARQUET_PATH)
        else:
            comments_df = pl.DataFrame({
                'id': pl.Series([], dtype=pl.Int64),
                'variant_key': pl.Series([], dtype=pl.Utf8),
                'sample_id': pl.Series([], dtype=pl.Utf8),
                'user_name': pl.Series([], dtype=pl.Utf8),
                'comment_text': pl.Series([], dtype=pl.Utf8),
                'timestamp': pl.Series([], dtype=pl.Utf8)
            })
        
        # Check DataFrame length properly
        new_id = comments_df['id'].max() + 1 if len(comments_df) > 0 else 1
        
        new_comment = pl.DataFrame({
            'id': [new_id],
            'variant_key': [variant_key],
            'sample_id': [sample_id],
            'user_name': [user_name],
            'comment_text': [comment_text],
            'timestamp': [datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
        })
        
        comments_df = pl.concat([comments_df, new_comment])
        comments_df.write_parquet(COMMENTS_PARQUET_PATH)
        
        # Invalidate cache
        db._invalidate_cache()
        
        return True
        
    except Exception as e:
        logger.error(f"Error adding comment: {e}")
        return False

def get_variant_comments(variant_key, sample_id):
    """Get all comments for a specific variant and sample"""
    try:
        if not os.path.exists(COMMENTS_PARQUET_PATH):
            return pd.DataFrame()
        
        comments_df = pl.read_parquet(COMMENTS_PARQUET_PATH)
        
        # Check DataFrame length properly
        if len(comments_df) == 0:
            return pd.DataFrame()
        
        variant_comments = comments_df.filter(
            (pl.col('variant_key') == variant_key) & 
            (pl.col('sample_id') == sample_id)
        ).sort('timestamp', descending=True)
        
        # Convert to pandas for easier handling in UI
        if len(variant_comments) > 0:
            return variant_comments.select(['user_name', 'comment_text', 'timestamp']).to_pandas()
        
        return pd.DataFrame()
        
    except Exception as e:
        logger.error(f"Error loading comments: {e}")
        return pd.DataFrame()

def update_variant_review_status(variant_key, sample_id, new_status):
    """Update the review status of a variant"""
    try:
        if not os.path.exists(VARIANTS_PARQUET_PATH):
            return False
        
        # Read the current data
        df = pl.read_parquet(VARIANTS_PARQUET_PATH)
        
        # Update the review status
        df = df.with_columns(
            pl.when(
                (pl.col('variant_key') == variant_key) & 
                (pl.col('SAMPLE') == sample_id)
            )
            .then(pl.lit(new_status))
            .otherwise(pl.col('review_status'))
            .alias('review_status')
        )
        
        # Write back to parquet
        df.write_parquet(VARIANTS_PARQUET_PATH)
        
        # Invalidate cache
        db._invalidate_cache()
        
        return True
        
    except Exception as e:
        logger.error(f"Error updating review status: {e}")
        return False

def get_variant_by_key(variant_key, sample_id):
    """Get a specific variant by its key and sample"""
    try:
        if not os.path.exists(VARIANTS_PARQUET_PATH):
            return None
        
        df = pl.read_parquet(VARIANTS_PARQUET_PATH)
        variant = df.filter(
            (pl.col('variant_key') == variant_key) & 
            (pl.col('SAMPLE') == sample_id)
        )
        
        if len(variant) > 0:
            # Ensure max_gnomad_af is available
            variant = db._ensure_max_gnomad_af(variant)
            return variant.to_pandas().iloc[0].to_dict()
        
        return None
        
    except Exception as e:
        logger.error(f"Error getting variant: {e}")
        return None

def export_variants_to_csv(df, filename=None):
    """Export variants DataFrame to CSV"""
    try:
        if filename is None:
            filename = f"variants_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        
        # Convert to pandas if it's Polars
        if isinstance(df, pl.DataFrame):
            df_pandas = df.to_pandas()
        else:
            df_pandas = df
        
        # Select relevant columns for export, including new frequency fields
        export_columns = [
            'CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT', 'VAF', 'gene', 
            'consequence', 'clinvar_sig', 'clinvar_id', 'clinvar_disease',
            'gnomad_af', 'max_gnomad_af', 'gnomad_af_afr', 'gnomad_af_amr', 
            'gnomad_af_asj', 'gnomad_af_eas', 'gnomad_af_fin', 'gnomad_af_nfe', 
            'gnomad_af_sas', 'ac_gnomad', 'nhomalt_gnomad', 'nhemalt_gnomad',
            'af_cgen', 'ac_cgen', 'an_cgen', 'cadd_score', 'review_status'
        ]
        
        available_columns = [col for col in export_columns if col in df_pandas.columns]
        export_df = df_pandas[available_columns].copy()
        
        # Save to CSV
        export_df.to_csv(filename, index=False)
        
        return filename
        
    except Exception as e:
        logger.error(f"Error exporting to CSV: {e}")
        return None

def get_database_info():
    """Get comprehensive database information"""
    try:
        stats = db.get_database_stats()
        
        # Add file information
        if os.path.exists(VARIANTS_PARQUET_PATH):
            file_size = os.path.getsize(VARIANTS_PARQUET_PATH) / (1024 * 1024)  # MB
            last_modified = datetime.fromtimestamp(os.path.getmtime(VARIANTS_PARQUET_PATH))
            
            stats.update({
                'file_size_mb': round(file_size, 2),
                'last_modified': last_modified.strftime("%Y-%m-%d %H:%M:%S"),
                'database_type': 'Parquet',
                'status': 'Connected'
            })
        else:
            stats.update({
                'file_size_mb': 0,
                'last_modified': 'N/A',
                'database_type': 'Parquet',
                'status': 'No data file found'
            })
        
        return stats
        
    except Exception as e:
        logger.error(f"Error getting database info: {e}")
        return {
            'total_variants': 0,
            'total_samples': 0,
            'reviewed_variants': 0,
            'pending_variants': 0,
            'clinvar_annotated': 0,
            'chromosomes': [],
            'file_size_mb': 0,
            'last_modified': 'Error',
            'database_type': 'Parquet',
            'status': 'Error'
        }