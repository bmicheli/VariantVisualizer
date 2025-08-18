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
            "gnomad_af": self._safe_float(annotations.get("gnomAD_genome_AF", "0")),
            "cadd_score": self._safe_float(annotations.get("CADD_phred", "")),
            "sift_score": self._safe_float(annotations.get("SIFT_score", "")),
            "polyphen_score": self._safe_float(
                annotations.get("Polyphen2_HVAR_score", "")
            ),
            "revel_score": self._safe_float(annotations.get("REVEL_score", "")),
            # Correct ClinVar fields
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
            "splice_ai": self._safe_float(annotations.get("spliceAI", "0")),
            "pli_score": self._safe_float(annotations.get("pLI", "")),
            "primateai_score": self._safe_float(annotations.get("primateAI", "")),
        }

        return result

    def _safe_float(self, value: str) -> Optional[float]:
        if not value or value == ".":
            return None
        try:
            return float(value.split(",")[0])
        except Exception:
            return None

    def _safe_int(self, value: str) -> Optional[int]:
        if not value or value == ".":
            return None
        try:
            return int(value.split(",")[0])
        except Exception:
            return None

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
        if ":" in aa_change:
            parts = aa_change.split(":")
            for part in parts:
                if part.startswith("p."):
                    return part
        return aa_change if aa_change.startswith("p.") else "p.?"

    def _parse_clinvar(self, clinvar: str) -> Optional[str]:
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
            return {"genotype": "./.", "depth": 0, "vaf": 0.0, "quality": 0}
        format_fields = format_str.split(":")
        sample_values = sample_str.split(":")
        while len(sample_values) < len(format_fields):
            sample_values.append(".")
        sample_data = dict(zip(format_fields, sample_values))
        return {
            "genotype": sample_data.get("GT", "./."),
            "depth": self._safe_int(sample_data.get("DP", "0")) or 0,
            "vaf": self._calculate_vaf(sample_data),
            "quality": self._safe_float(sample_data.get("GQ", "0")) or 0,
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
                    if sample_info["genotype"] in ["0/0", "0|0", "./."]:
                        continue
                    record = {
                        "CHROM": chrom.replace("chr", ""),
                        "POS": int(pos),
                        "ID": var_id if var_id != "." else "",
                        "REF": ref,
                        "ALT": alt,
                        "QUAL": self._safe_float(qual) or 0,
                        "FILTER": filter_val,
                        "SAMPLE": sample,
                        "GT": sample_info["genotype"],
                        "DP": sample_info["depth"],
                        "VAF": sample_info["vaf"],
                        "GQ": sample_info["quality"],
                        "AD": sample_info["allelic_depth"],
                        "variant_key": variant_key,
                        "gene": str(annotations.get("gene", "UNKNOWN")),
                        "consequence": str(annotations.get("consequence", "variant")),
                        "aa_change": str(annotations.get("aa_change", "p.?")),
                        "clinvar_sig": annotations.get("clinvar_sig") or "",
                        "clinvar_id": annotations.get("clinvar_id") or "",
                        "clinvar_disease": annotations.get("clinvar_disease") or "",
                        "clinvar_review_status": annotations.get("clinvar_review_status") or "",
                        "upload_date": datetime.now().isoformat(),
                        "review_status": "Pending",
                        **{k: v for k, v in annotations.items() if k not in [
                            "gene","consequence","aa_change","clinvar_sig","clinvar_id","clinvar_disease","clinvar_review_status"
                        ]},
                    }
                    variant_records.append(record)
        if not variant_records:
            return pl.DataFrame()
        df = pl.DataFrame(variant_records)
        df = df.with_columns([
            pl.col("CHROM").cast(pl.Utf8),
            pl.col("POS").cast(pl.UInt32),
            pl.col("QUAL").cast(pl.Float32),
            pl.col("DP").cast(pl.UInt16),
            pl.col("VAF").cast(pl.Float32),
            pl.col("GQ").cast(pl.Float32),
            pl.col("af").cast(pl.Float32),
            pl.col("gnomad_af").cast(pl.Float32),
            pl.col("cadd_score").cast(pl.Float32),
            pl.col("sift_score").cast(pl.Float32),
            pl.col("polyphen_score").cast(pl.Float32),
            pl.col("splice_ai").cast(pl.Float32),
            pl.col("gene").cast(pl.Utf8),
            pl.col("consequence").cast(pl.Utf8),
            pl.col("aa_change").cast(pl.Utf8),
            pl.col("clinvar_sig").cast(pl.Utf8),
            pl.col("clinvar_id").cast(pl.Utf8),
            pl.col("clinvar_disease").cast(pl.Utf8),
            pl.col("clinvar_review_status").cast(pl.Utf8),
        ])
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
                    chunk_lines = []
            if chunk_lines:
                chunk_df = self.process_vcf_chunk(chunk_lines, samples)
                if not chunk_df.is_empty():
                    chunk_dfs.append(chunk_df)
                    total_variants += len(chunk_df)
        if not chunk_dfs:
            logger.error("No variant data found in VCF file")
            return ""
        final_df = pl.concat(chunk_dfs)
        final_df = final_df.sort(["CHROM", "POS"])
        final_df.write_parquet(
            output_file,
            compression="snappy",
            use_pyarrow=True,
            row_group_size=50000,
            statistics=True,
        )
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
	
	args = parser.parse_args()
	
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