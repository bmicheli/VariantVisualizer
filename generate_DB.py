#!/usr/bin/env python3
import polars as pl
import pandas as pd
import numpy as np
import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
import argparse
import sys

def check_dependencies():
	"""Check if required dependencies are available"""
	missing_deps = []
	
	try:
		import pyarrow
		import pyarrow.parquet
		logger.info(f"PyArrow version: {pyarrow.__version__}")
	except ImportError:
		missing_deps.append("pyarrow")
		logger.warning("PyArrow not found - will use native parquet writer")
	
	if missing_deps:
		logger.info("To install missing dependencies, run:")
		logger.info(f"pip install {' '.join(missing_deps)}")
	
	return len(missing_deps) == 0

# Configure logging
logging.basicConfig(
	level=logging.INFO,
	format="%(asctime)s - %(levelname)s - %(message)s",
	handlers=[
		logging.FileHandler("vcf_conversion.log"),
		logging.StreamHandler(sys.stdout),
	],
)
logger = logging.getLogger(__name__)


class VCFToParquetConverter:
	def __init__(self, output_dir: str = "data", chunk_size: int = 10000):
		self.output_dir = Path(output_dir)
		self.chunk_size = chunk_size
		self.output_dir.mkdir(exist_ok=True)

		# ClinVar significance mappings
		self.clinvar_mapping = {
			"pathogenic": "Pathogenic",
			"likely_pathogenic": "Likely pathogenic",
			"uncertain_significance": "VUS",
			"likely_benign": "Likely benign",
			"benign": "Benign",
			"conflicting_interpretations": "Conflicting",
			"conflicting_interpretations_of_pathogenicity": "Conflicting",
			"other": "Other",
			"not_provided": "Not provided",
			"drug_response": "Drug response",
			"association": "Association",
			"protective": "Protective",
			"affects": "Affects",
			"risk_factor": "Risk factor",
		}

		# Consequence severity mapping
		self.consequence_severity = {
			"transcript_ablation": "high",
			"splice_acceptor_variant": "high",
			"splice_donor_variant": "high",
			"stop_gained": "high",
			"frameshift_variant": "high",
			"stop_lost": "high",
			"start_lost": "high",
			"transcript_amplification": "high",
			"inframe_insertion": "moderate",
			"inframe_deletion": "moderate",
			"missense_variant": "moderate",
			"protein_altering_variant": "moderate",
			"splice_region_variant": "low",
			"incomplete_terminal_codon_variant": "low",
			"start_retained_variant": "low",
			"stop_retained_variant": "low",
			"synonymous_variant": "low",
			"coding_sequence_variant": "modifier",
			"mature_miRNA_variant": "modifier",
			"5_prime_UTR_variant": "modifier",
			"3_prime_UTR_variant": "modifier",
			"non_coding_transcript_exon_variant": "modifier",
			"intron_variant": "modifier",
			"NMD_transcript_variant": "modifier",
			"non_coding_transcript_variant": "modifier",
			"upstream_gene_variant": "modifier",
			"downstream_gene_variant": "modifier",
			"TFBS_ablation": "modifier",
			"TFBS_amplification": "modifier",
			"TF_binding_site_variant": "modifier",
			"regulatory_region_ablation": "modifier",
			"regulatory_region_amplification": "modifier",
			"feature_elongation": "modifier",
			"regulatory_region_variant": "modifier",
			"feature_truncation": "modifier",
			"intergenic_variant": "modifier",
		}

	def parse_vcf_header(self, vcf_file: str) -> Tuple[List[str], Dict[str, Any]]:
		"""Parse VCF header to extract sample names and metadata"""
		samples = []
		info_fields = {}
		format_fields = {}

		with open(vcf_file, "r") as f:
			for line in f:
				if line.startswith("##INFO"):
					match = re.search(r"ID=([^,]+)", line)
					if match:
						field_id = match.group(1)
						info_fields[field_id] = line.strip()
				elif line.startswith("##FORMAT"):
					match = re.search(r"ID=([^,]+)", line)
					if match:
						field_id = match.group(1)
						format_fields[field_id] = line.strip()
				elif line.startswith("#CHROM"):
					headers = line.strip().split("\t")
					if len(headers) > 9:
						samples = headers[9:]
					break
				elif not line.startswith("#"):
					break

		logger.info(f"Found {len(samples)} samples in VCF file")
		return samples, {"info": info_fields, "format": format_fields}

	def extract_info_annotations(self, info_str: str) -> Dict[str, Any]:
		"""Extract and parse annotations from INFO field"""
		annotations = {}

		if not info_str or info_str == ".":
			return annotations

		for item in info_str.split(";"):
			if "=" in item:
				key, value = item.split("=", 1)
				annotations[key] = value
			else:
				annotations[item] = True

		# Parse gnomAD population-specific frequencies
		gnomad_af_afr = self._safe_float(annotations.get("AF_gnomad_afr", "0"))
		gnomad_af_amr = self._safe_float(annotations.get("AF_gnomad_amr", "0"))
		gnomad_af_eas = self._safe_float(annotations.get("AF_gnomad_eas", "0"))
		gnomad_af_nfe = self._safe_float(annotations.get("AF_gnomad_nfe", "0"))
		gnomad_af_sas = self._safe_float(annotations.get("AF_gnomad_sas", "0"))
		gnomad_af_asj = self._safe_float(annotations.get("AF_gnomad_asj", "0"))
		gnomad_af_fin = self._safe_float(annotations.get("AF_gnomad_fin", "0"))

		# Calculate max gnomAD AF
		gnomad_afs = [gnomad_af_afr, gnomad_af_amr, gnomad_af_eas, gnomad_af_nfe, gnomad_af_sas,gnomad_af_asj,gnomad_af_fin]
		max_gnomad_af = max([af for af in gnomad_afs if af > 0], default=0.0)

		result = {
			"af": self._safe_float(annotations.get("AF", "0")),
			"ac": self._safe_int(annotations.get("AC", "0")),
			"an": self._safe_int(annotations.get("AN", "0")),
			"qual": self._safe_float(annotations.get("AQ", "0")),
			"gene": annotations.get("Gene.HGNC", "UNKNOWN"),
			"consequence": self._parse_consequence(
				annotations.get("ExonicFunc.HGNC", "variant")
			),
			"aa_change": self._parse_aa_change(annotations.get("AAChange.HGNC", "")),
			
			# Original gnomAD fields
			"gnomad_af": self._safe_float(annotations.get("AF_gnomad", "0")),
			
			# New gnomAD population-specific fields
			"gnomad_af_afr": gnomad_af_afr,
			"gnomad_af_amr": gnomad_af_amr,
			"gnomad_af_asj": gnomad_af_asj,
			"gnomad_af_eas": gnomad_af_eas,
			"gnomad_af_fin": gnomad_af_fin,
			"gnomad_af_nfe": gnomad_af_nfe,
			"gnomad_af_sas": gnomad_af_sas,
			"max_gnomad_af": max_gnomad_af,  # This will be used for display
			
			# gnomAD counts
			"ac_gnomad": self._safe_int(annotations.get("AC_gnomad", "0")),
			"nhomalt_gnomad": self._safe_float(annotations.get("nhomalt_gnomad", "0")),
			"nhemalt_gnomad": self._safe_float(annotations.get("nhemalt_gnomad", "0")),
			
			# CGEN frequencies
			"ac_cgen": self._safe_int(annotations.get("AC_CGEN", "0")),
			"af_cgen": self._safe_float(annotations.get("AF_CGEN", "0")),
			"an_cgen": self._safe_int(annotations.get("AN_CGEN", "0")),
			
			# Prediction scores
			"cadd_score": self._safe_float(annotations.get("CADD_phred", "")),
			"sift_score": self._safe_float(annotations.get("SIFT_score", "")),
			"polyphen_score": self._safe_float(
				annotations.get("Polyphen2_HVAR_score", "")
			),
			"revel_score": self._safe_float(annotations.get("REVEL_score", "")),
			
			# ClinVar fields
			"clinvar_sig": self._parse_clinvar(
				annotations.get("CLNSIG", "") or annotations.get("CLNSIGCONF", "")
			),
			"clinvar_id": annotations.get("CLNVID")
			or annotations.get("ClinVar_ID")
			or "",
			"clinvar_disease": annotations.get("CLNDBN")
			or annotations.get("ClinVar_CLNDBN")
			or "",
			"clinvar_review_status": annotations.get("CLNREVSTAT")
			or annotations.get("ClinVar_CLNREVSTAT")
			or "",
			
			# Additional scores
			"splice_ai": self._safe_float(annotations.get("spliceAI", "0")),
			"pli_score": self._safe_float(annotations.get("pLI", "")),
			"primateai_score": self._safe_float(annotations.get("primateAI", "")),
		}

		return result

	def _safe_float(self, value: str) -> float:
		"""Safely convert string to float, return 0.0 for invalid values"""
		if not value or value == "." or value == "":
			return 0.0
		try:
			# Handle space-separated values by taking the first one
			if " " in str(value):
				value = str(value).split()[0]
			return float(value.split(",")[0])
		except (ValueError, TypeError):
			return 0.0

	def _safe_int(self, value: str) -> int:
		"""Safely convert string to int, return 0 for invalid values"""
		if not value or value == "." or value == "":
			return 0
		try:
			# Handle space-separated values by taking the first one
			if " " in str(value):
				value = str(value).split()[0]
			return int(float(value.split(",")[0]))
		except (ValueError, TypeError):
			return 0

	def _parse_consequence(self, consequence: str) -> str:
		if not consequence or consequence == ".":
			return "variant"
		consequence_map = {
			"nonsynonymous_SNV": "missense_variant",
			"synonymous_SNV": "synonymous_variant",
			"stopgain": "stop_gained",
			"stoploss": "stop_lost",
			"frameshift_deletion": "frameshift_variant",
			"frameshift_insertion": "frameshift_variant",
		}
		return consequence_map.get(consequence, consequence)

	def _parse_aa_change(self, aa_change: str) -> str:
		if not aa_change or aa_change == ".":
			return "p.?"
		# Keep only entries that contain a protein change (p.)
		entries = [part for part in aa_change.split(",") if "p." in part]
		return ",".join(entries) if entries else "p.?"

	def _parse_clinvar(self, clinvar: str) -> str:
		"""Parse ClinVar significance, return empty string for invalid values"""
		if not clinvar or clinvar == ".":
			return ""
		clinvar_values = re.split(r"[,;|/]", clinvar.lower())
		priority_order = [
			"pathogenic",
			"likely_pathogenic",
			"conflicting_interpretations",
			"uncertain_significance",
			"likely_benign",
			"benign",
		]
		for priority_val in priority_order:
			for val in clinvar_values:
				if priority_val in val:
					return self.clinvar_mapping.get(priority_val, "Other")
		for val in clinvar_values:
			for key, mapped_value in self.clinvar_mapping.items():
				if key in val:
					return mapped_value
		return "Other"

	def parse_sample_data(self, format_str: str, sample_str: str) -> Dict[str, Any]:
		if not format_str or not sample_str or sample_str == ".":
			return {"genotype": "./.", "depth": 0, "vaf": 0.0, "quality": 0.0, "allelic_depth": "0,0"}
		
		format_fields = format_str.split(":")
		sample_values = sample_str.split(":")
		
		# Pad sample values with "." if needed
		while len(sample_values) < len(format_fields):
			sample_values.append(".")
			
		sample_data = dict(zip(format_fields, sample_values))
		
		return {
			"genotype": sample_data.get("GT", "./."),
			"depth": self._safe_int(sample_data.get("DP", "0")),
			"vaf": self._calculate_vaf(sample_data),
			"quality": self._safe_float(sample_data.get("GQ", "0")),
			"allelic_depth": sample_data.get("AD", "0,0"),
		}

	def _calculate_vaf(self, sample_data: Dict[str, str]) -> float:
		if "VAF" in sample_data and sample_data["VAF"] != ".":
			try:
				return float(sample_data["VAF"])
			except Exception:
				pass
		ad = sample_data.get("AD", "0,0")
		if ad != "." and "," in ad:
			try:
				depths = [int(x) for x in ad.split(",")]
				if len(depths) >= 2 and sum(depths) > 0:
					return depths[1] / sum(depths)
			except Exception:
				pass
		return 0.0

	def create_empty_dataframe_with_schema(self) -> pl.DataFrame:
		"""Create an empty DataFrame with the correct schema"""
		schema = {
			"CHROM": pl.Utf8,
			"POS": pl.UInt32,
			"ID": pl.Utf8,
			"REF": pl.Utf8,
			"ALT": pl.Utf8,
			"QUAL": pl.Float32,
			"FILTER": pl.Utf8,
			"SAMPLE": pl.Utf8,
			"GT": pl.Utf8,
			"DP": pl.UInt16,
			"VAF": pl.Float32,
			"GQ": pl.Float32,
			"AD": pl.Utf8,
			"variant_key": pl.Utf8,
			"gene": pl.Utf8,
			"consequence": pl.Utf8,
			"aa_change": pl.Utf8,
			"clinvar_sig": pl.Utf8,
			"clinvar_id": pl.Utf8,
			"clinvar_disease": pl.Utf8,
			"clinvar_review_status": pl.Utf8,
			"upload_date": pl.Utf8,
			"review_status": pl.Utf8,
			"af": pl.Float32,
			"ac": pl.Int32,
			"an": pl.Int32,
			"qual": pl.Float32,
			"gnomad_af": pl.Float32,
			"gnomad_af_afr": pl.Float32,
			"gnomad_af_amr": pl.Float32,
			"gnomad_af_asj": pl.Float32,
			"gnomad_af_eas": pl.Float32,
			"gnomad_af_fin": pl.Float32,
			"gnomad_af_nfe": pl.Float32,
			"gnomad_af_sas": pl.Float32,
			"max_gnomad_af": pl.Float32,
			"ac_gnomad": pl.Int32,
			"nhomalt_gnomad": pl.Float32,
			"nhemalt_gnomad": pl.Float32,
			"ac_cgen": pl.Int32,
			"af_cgen": pl.Float32,
			"an_cgen": pl.Int32,
			"cadd_score": pl.Float32,
			"sift_score": pl.Float32,
			"polyphen_score": pl.Float32,
			"revel_score": pl.Float32,
			"splice_ai": pl.Float32,
			"pli_score": pl.Float32,
			"primateai_score": pl.Float32,
		}
		return pl.DataFrame(schema=schema)

	def process_vcf_chunk(self, chunk_lines: List[str], samples: List[str]) -> pl.DataFrame:
		variant_records = []
		
		for line in chunk_lines:
			if line.startswith("#") or not line.strip():
				continue
			
			fields = line.strip().split("\t")
			if len(fields) < 8:
				continue
				
			chrom, pos, var_id, ref, alt, qual, filter_val, info = fields[:8]
			format_str = fields[8] if len(fields) > 8 else ""
			sample_data = fields[9:] if len(fields) > 9 else []
			
			annotations = self.extract_info_annotations(info)
			variant_key = f"{chrom}:{pos}:{ref}:{alt}"
			
			for i, sample in enumerate(samples):
				if i < len(sample_data):
					sample_info = self.parse_sample_data(format_str, sample_data[i])
					
					# Skip reference calls
					if sample_info["genotype"] in ["0/0", "0|0", "./."]:
						continue
					
					# Create record with consistent types
					record = {
						"CHROM": str(chrom.replace("chr", "")),
						"POS": int(pos),
						"ID": str(var_id) if var_id != "." else "",
						"REF": str(ref),
						"ALT": str(alt),
						"QUAL": self._safe_float(qual),
						"FILTER": str(filter_val),
						"SAMPLE": str(sample),
						"GT": str(sample_info["genotype"]),
						"DP": int(sample_info["depth"]),
						"VAF": float(sample_info["vaf"]),
						"GQ": float(sample_info["quality"]),
						"AD": str(sample_info["allelic_depth"]),
						"variant_key": str(variant_key),
						"gene": str(annotations.get("gene", "UNKNOWN")),
						"consequence": str(annotations.get("consequence", "variant")),
						"aa_change": str(annotations.get("aa_change", "p.?")),
						"clinvar_sig": str(annotations.get("clinvar_sig", "")),
						"clinvar_id": str(annotations.get("clinvar_id", "")),
						"clinvar_disease": str(annotations.get("clinvar_disease", "")),
						"clinvar_review_status": str(annotations.get("clinvar_review_status", "")),
						"upload_date": str(datetime.now().isoformat()),
						"review_status": "Pending",
						"af": float(annotations.get("af", 0.0)),
						"ac": int(annotations.get("ac", 0)),
						"an": int(annotations.get("an", 0)),
						"qual": float(annotations.get("qual", 0.0)),
						"gnomad_af": float(annotations.get("gnomad_af", 0.0)),
						"gnomad_af_afr": float(annotations.get("gnomad_af_afr", 0.0)),
						"gnomad_af_amr": float(annotations.get("gnomad_af_amr", 0.0)),
						"gnomad_af_asj": float(annotations.get("gnomad_af_asj", 0.0)),
						"gnomad_af_eas": float(annotations.get("gnomad_af_eas", 0.0)),
						"gnomad_af_fin": float(annotations.get("gnomad_af_fin", 0.0)),
						"gnomad_af_nfe": float(annotations.get("gnomad_af_nfe", 0.0)),
						"gnomad_af_sas": float(annotations.get("gnomad_af_sas", 0.0)),
						"max_gnomad_af": float(annotations.get("max_gnomad_af", 0.0)),
						"ac_gnomad": int(annotations.get("ac_gnomad", 0)),
						"nhomalt_gnomad": float(annotations.get("nhomalt_gnomad", 0.0)),
						"nhemalt_gnomad": float(annotations.get("nhemalt_gnomad", 0.0)),
						"ac_cgen": int(annotations.get("ac_cgen", 0)),
						"af_cgen": float(annotations.get("af_cgen", 0.0)),
						"an_cgen": int(annotations.get("an_cgen", 0)),
						"cadd_score": float(annotations.get("cadd_score", 0.0)),
						"sift_score": float(annotations.get("sift_score", 0.0)),
						"polyphen_score": float(annotations.get("polyphen_score", 0.0)),
						"revel_score": float(annotations.get("revel_score", 0.0)),
						"splice_ai": float(annotations.get("splice_ai", 0.0)),
						"pli_score": float(annotations.get("pli_score", 0.0)),
						"primateai_score": float(annotations.get("primateai_score", 0.0)),
					}
					variant_records.append(record)
		
		if not variant_records:
			return self.create_empty_dataframe_with_schema()
		
		# Create DataFrame with explicit schema
		df = pl.DataFrame(variant_records, schema={
			"CHROM": pl.Utf8,
			"POS": pl.UInt32,
			"ID": pl.Utf8,
			"REF": pl.Utf8,
			"ALT": pl.Utf8,
			"QUAL": pl.Float32,
			"FILTER": pl.Utf8,
			"SAMPLE": pl.Utf8,
			"GT": pl.Utf8,
			"DP": pl.UInt16,
			"VAF": pl.Float32,
			"GQ": pl.Float32,
			"AD": pl.Utf8,
			"variant_key": pl.Utf8,
			"gene": pl.Utf8,
			"consequence": pl.Utf8,
			"aa_change": pl.Utf8,
			"clinvar_sig": pl.Utf8,
			"clinvar_id": pl.Utf8,
			"clinvar_disease": pl.Utf8,
			"clinvar_review_status": pl.Utf8,
			"upload_date": pl.Utf8,
			"review_status": pl.Utf8,
			"af": pl.Float32,
			"ac": pl.Int32,
			"an": pl.Int32,
			"qual": pl.Float32,
			"gnomad_af": pl.Float32,
			"gnomad_af_afr": pl.Float32,
			"gnomad_af_amr": pl.Float32,
			"gnomad_af_asj": pl.Float32,
			"gnomad_af_eas": pl.Float32,
			"gnomad_af_fin": pl.Float32,
			"gnomad_af_nfe": pl.Float32,
			"gnomad_af_sas": pl.Float32,
			"max_gnomad_af": pl.Float32,
			"ac_gnomad": pl.Int32,
			"nhomalt_gnomad": pl.Float32,
			"nhemalt_gnomad": pl.Float32,
			"ac_cgen": pl.Int32,
			"af_cgen": pl.Float32,
			"an_cgen": pl.Int32,
			"cadd_score": pl.Float32,
			"sift_score": pl.Float32,
			"polyphen_score": pl.Float32,
			"revel_score": pl.Float32,
			"splice_ai": pl.Float32,
			"pli_score": pl.Float32,
			"primateai_score": pl.Float32,
		})
		
		return df

	def convert_vcf_to_parquet(self, vcf_file: str, output_name: str = "variants") -> str:
		logger.info(f"Starting conversion of {vcf_file}")
		samples, metadata = self.parse_vcf_header(vcf_file)
		output_file = self.output_dir / f"{output_name}.parquet"
		
		chunk_dfs = []
		total_variants = 0
		
		with open(vcf_file, "r") as f:
			chunk_lines = []
			line_count = 0
			
			for line in f:
				if line.startswith("#"):
					continue
					
				chunk_lines.append(line)
				line_count += 1
				
				if len(chunk_lines) >= self.chunk_size:
					chunk_df = self.process_vcf_chunk(chunk_lines, samples)
					if not chunk_df.is_empty():
						chunk_dfs.append(chunk_df)
						total_variants += len(chunk_df)
						logger.info(f"Processed {line_count} lines, {total_variants} variants so far")
					chunk_lines = []
			
			# Process remaining lines
			if chunk_lines:
				chunk_df = self.process_vcf_chunk(chunk_lines, samples)
				if not chunk_df.is_empty():
					chunk_dfs.append(chunk_df)
					total_variants += len(chunk_df)
		
		if not chunk_dfs:
			logger.error("No variant data found in VCF file")
			return ""
		
		# Concatenate all chunks
		logger.info("Concatenating all chunks...")
		final_df = pl.concat(chunk_dfs)
		final_df = final_df.sort(["CHROM", "POS"])
		
		# Write to parquet
		logger.info("Writing to parquet file...")
		try:
			# Try with PyArrow first
			final_df.write_parquet(
				output_file,
				compression="snappy",
				use_pyarrow=True,
				row_group_size=50000,
				statistics=True,
			)
		except (ModuleNotFoundError, ImportError) as e:
			logger.warning(f"PyArrow not available ({e}), falling back to native parquet writer")
			try:
				final_df.write_parquet(
					output_file,
					compression="snappy",
					use_pyarrow=False,
				)
			except Exception as e:
				logger.error(f"Failed to write parquet file: {e}")
				# Fallback to CSV if parquet fails
				csv_output = output_file.with_suffix('.csv')
				logger.info(f"Falling back to CSV format: {csv_output}")
				final_df.write_csv(csv_output)
				return str(csv_output)
		
		# Create metadata
		stats = {
			"total_variants": len(final_df),
			"total_samples": len(samples),
			"chromosomes": final_df["CHROM"].unique().to_list(),
			"clinvar_annotated": len(final_df.filter(pl.col("clinvar_sig") != "")),
			"date_created": datetime.now().isoformat(),
			"source_file": vcf_file,
		}
		
		metadata_file = self.output_dir / f"{output_name}_metadata.json"
		import json
		with open(metadata_file, "w") as f:
			json.dump(stats, f, indent=2)
		
		logger.info(f"Conversion complete: {total_variants} variants saved to {output_file}")
		return str(output_file)

	def create_sample_index(self, parquet_file: str):
		logger.info("Creating sample index...")
		df = pl.read_parquet(parquet_file)
		
		sample_stats = (
			df.group_by("SAMPLE")
			.agg([
				pl.count("variant_key").alias("variant_count"),
				pl.col("CHROM").n_unique().alias("chromosomes_with_variants"),
				pl.col("consequence").value_counts().alias("consequences"),
				pl.col("clinvar_sig").filter(pl.col("clinvar_sig") != "").len().alias("clinvar_annotated"),
			])
		)
		
		index_file = self.output_dir / "sample_index.parquet"
		sample_stats.write_parquet(index_file)
		logger.info(f"Sample index created: {index_file}")

def main():
	"""Main execution function"""
	parser = argparse.ArgumentParser(description="Convert VCF files to optimized Parquet format (ClinVar focused)")
	parser.add_argument("vcf_file", help="Input VCF file path")
	parser.add_argument("-o", "--output-dir", default="data", help="Output directory")
	parser.add_argument("-n", "--name", default="variants", help="Output file name prefix")
	parser.add_argument("-c", "--chunk-size", type=int, default=10000, help="Chunk size for processing")
	parser.add_argument("--create-index", action="store_true", help="Create sample index")
	parser.add_argument("--skip-deps-check", action="store_true", help="Skip dependency check")
	
	args = parser.parse_args()
	
	# Check dependencies unless explicitly skipped
	if not args.skip_deps_check:
		check_dependencies()
	
	if not os.path.exists(args.vcf_file):
		logger.error(f"VCF file not found: {args.vcf_file}")
		sys.exit(1)
	
	# Initialize converter
	converter = VCFToParquetConverter(args.output_dir, args.chunk_size)
	
	# Convert VCF to Parquet
	output_file = converter.convert_vcf_to_parquet(args.vcf_file, args.name)
	
	if output_file and args.create_index:
		converter.create_sample_index(output_file)
	
	logger.info("Conversion pipeline completed successfully!")

if __name__ == "__main__":
	main()

# Example usage:
# python vcf_to_parquet.py merged_variants.vcf -o data -n variants --create-index
# python vcf_to_parquet.py large_cohort.vcf -c 20000 --create-index