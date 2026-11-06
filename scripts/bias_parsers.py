"""
Parsers for the BIAS-2015 supplementary files used by parse_bias_data.py:
the eRepo VCF, InterVar predictions, BIAS-2015 predictions, and the
benchmark ID alignment map. Also HGVS and VCF coordinate helpers.
"""

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

# Increase CSV field size limit for large annotation fields
csv.field_size_limit(sys.maxsize)


def load_id_alignment_map(original_dataset_path: Path):
    """
    Creates a map from variation/biological key to the original benchmark's UUID.
    This handles cases where the source files have different UUIDs for the same variant.
    """
    if not original_dataset_path.exists():
        print(f"Warning: Original dataset not found at {original_dataset_path}. Alignment skipping.")
        return {}, set()
    
    alignment_map = {} # variation -> original_id
    original_ids = set()
    with open(original_dataset_path, 'r', encoding='utf-8') as f:
        for line in f:
            data = json.loads(line)
            oid = data.get('id')
            variation = str(data.get('variation', '')).lower()
            variation_orig = str(data.get('variation_original', '')).lower()
            
            if oid:
                original_ids.add(oid)
                # Map by UUID primarily
                alignment_map[oid] = oid
                # Map by variation strings for rescue remapping
                if variation: alignment_map[variation] = oid
                if variation_orig: alignment_map[variation_orig] = oid
    
    return alignment_map, original_ids


# HG19/GRCh37 specific NC_ accession to chromosome mapping
HG19_NC_TO_CHROM = {
    "NC_000001.10": "chr1", "NC_000001.9": "chr1",
    "NC_000002.11": "chr2", "NC_000002.10": "chr2",
    "NC_000003.11": "chr3", "NC_000003.10": "chr3",
    "NC_000004.11": "chr4", "NC_000004.10": "chr4",
    "NC_000005.9": "chr5", "NC_000005.8": "chr5",
    "NC_000006.11": "chr6", "NC_000006.10": "chr6",
    "NC_000007.13": "chr7", "NC_000007.12": "chr7",
    "NC_000008.10": "chr8", "NC_000008.9": "chr8",
    "NC_000009.11": "chr9", "NC_000009.10": "chr9",
    "NC_000010.10": "chr10", "NC_000010.9": "chr10",
    "NC_000011.9": "chr11", "NC_000011.8": "chr11",
    "NC_000012.11": "chr12", "NC_000012.10": "chr12",
    "NC_000013.10": "chr13", "NC_000013.9": "chr13",
    "NC_000014.8": "chr14", "NC_000014.7": "chr14",
    "NC_000015.9": "chr15", "NC_000015.8": "chr15",
    "NC_000016.9": "chr16", "NC_000016.8": "chr16",
    "NC_000017.10": "chr17", "NC_000017.9": "chr17",
    "NC_000018.9": "chr18", "NC_000018.8": "chr18",
    "NC_000019.9": "chr19", "NC_000019.8": "chr19",
    "NC_000020.10": "chr20", "NC_000020.9": "chr20",
    "NC_000021.8": "chr21", "NC_000021.7": "chr21",
    "NC_000022.10": "chr22", "NC_000022.9": "chr22",
    "NC_000023.10": "chrX", "NC_000023.9": "chrX",
    "NC_000024.9": "chrY", "NC_000024.8": "chrY",
}

NC_TO_CHROM = dict(HG19_NC_TO_CHROM)
for chrom_num in range(1, 23):
    for version in range(7, 15):
        key = f"NC_{chrom_num:06d}.{version}"
        if key not in NC_TO_CHROM:
            NC_TO_CHROM[key] = f"chr{chrom_num}"
for version in range(7, 15):
    key = f"NC_000023.{version}"
    if key not in NC_TO_CHROM:
        NC_TO_CHROM[key] = "chrX"
for version in range(7, 15):
    key = f"NC_000024.{version}"
    if key not in NC_TO_CHROM:
        NC_TO_CHROM[key] = "chrY"

def chrom_num_to_name(chrom_num: int) -> str:
    if chrom_num == 23:
        return "chrX"
    elif chrom_num == 24:
        return "chrY"
    else:
        return f"chr{chrom_num}"

def vcf_to_annovar_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return f'{chrom}:{pos}:{ref}:{alt}'
    prefix_len = 0
    for i in range(min(len(ref), len(alt))):
        if ref[i] == alt[i]:
            prefix_len += 1
        else:
            break
    new_ref = ref[prefix_len:] if prefix_len < len(ref) else ""
    new_alt = alt[prefix_len:] if prefix_len < len(alt) else ""
    new_pos = pos + prefix_len
    if new_alt == "" or (len(new_ref) > 0 and len(new_alt) == 0):
        return f'{chrom}:{new_pos}:{new_ref}:-'
    if new_ref == "" or (len(new_ref) == 0 and len(new_alt) > 0):
        return f'{chrom}:{new_pos - 1}:-:{new_alt}'
    return f'{chrom}:{new_pos}:{new_ref}:{new_alt}'

def vcf_to_annovar_keys(chrom: str, pos: int, ref: str, alt: str) -> list[str]:
    keys = [vcf_to_annovar_key(chrom, pos, ref, alt)]
    if len(ref) == len(alt) and len(ref) > 1:
        for i in range(len(ref)):
            if ref[i] != alt[i]:
                snv_key = f'{chrom}:{pos + i}:{ref[i]}:{alt[i]}'
                keys.append(snv_key)
    return keys

def parse_hg19_coord(hgvs_list: str) -> tuple[str, int, str, str, str] | None:
    for expr in hgvs_list.split(", "):
        expr = expr.strip()
        nc_match = re.match(r'(NC_\d+\.\d+):g\.(.+)', expr)
        if not nc_match:
            continue
        accession = nc_match.group(1)
        gen_coords = nc_match.group(2)
        if accession not in HG19_NC_TO_CHROM:
            continue
        chrom = HG19_NC_TO_CHROM[accession]
        snv_match = re.match(r'(\d+)([ACGT]+)>([ACGT]+)$', gen_coords)
        if snv_match:
            pos = int(snv_match.group(1))
            ref = snv_match.group(2)
            alt = snv_match.group(3)
            return (chrom, pos, ref, alt, expr)
        del_match = re.match(r'(\d+)(?:_(\d+))?del([ACGT]*)$', gen_coords)
        if del_match:
            pos = int(del_match.group(1))
            ref = del_match.group(3) if del_match.group(3) else "-"
            return (chrom, pos, ref, "-", expr)
        ins_match = re.match(r'(\d+)_\d+ins([ACGT]+)$', gen_coords)
        if ins_match:
            pos = int(ins_match.group(1))
            alt = ins_match.group(2)
            return (chrom, pos, "-", alt, expr)
        delins_match = re.match(r'(\d+)(?:_\d+)?delins([ACGT]+)$', gen_coords)
        if delins_match:
            pos = int(delins_match.group(1))
            alt = delins_match.group(2)
            return (chrom, pos, "delins", alt, expr)
        dup_match = re.match(r'(\d+)(?:_\d+)?dup([ACGT]*)$', gen_coords)
        if dup_match:
            pos = int(dup_match.group(1))
            dup_seq = dup_match.group(2) if dup_match.group(2) else "dup"
            return (chrom, pos, "-", dup_seq, expr)
    return None

def parse_intervar_classification(intervar_field: str) -> tuple[str, list[str]]:
    match = re.search(r'InterVar:\s*(\w+(?:\s+\w+)?)', intervar_field)
    if not match:
        return ("", [])
    classification = match.group(1).strip()
    classification_lower = classification.lower()
    if "pathogenic" in classification_lower and "likely" not in classification_lower:
        classification = "Pathogenic"
    elif "likely pathogenic" in classification_lower:
        classification = "Likely Pathogenic"
    elif "uncertain" in classification_lower:
        classification = "Uncertain Significance"
    elif "likely benign" in classification_lower:
        classification = "Likely Benign"
    elif "benign" in classification_lower:
        classification = "Benign"
    criteria_met = []
    if re.search(r'PVS1=1', intervar_field): criteria_met.append("PVS1")
    ps_match = re.search(r'PS=\[([^\]]+)\]', intervar_field)
    if ps_match:
        ps_vals = [int(x.strip()) for x in ps_match.group(1).split(',')]
        for i, val in enumerate(ps_vals):
            if val == 1: criteria_met.append(f"PS{i+1}")
    pm_match = re.search(r'PM=\[([^\]]+)\]', intervar_field)
    if pm_match:
        pm_vals = [int(x.strip()) for x in pm_match.group(1).split(',')]
        for i, val in enumerate(pm_vals):
            if val == 1: criteria_met.append(f"PM{i+1}")
    pp_match = re.search(r'PP=\[([^\]]+)\]', intervar_field)
    if pp_match:
        pp_vals = [int(x.strip()) for x in pp_match.group(1).split(',')]
        for i, val in enumerate(pp_vals):
            if val == 1: criteria_met.append(f"PP{i+1}")
    if re.search(r'BA1=1', intervar_field): criteria_met.append("BA1")
    bs_match = re.search(r'BS=\[([^\]]+)\]', intervar_field)
    if bs_match:
        bs_vals = [int(x.strip()) for x in bs_match.group(1).split(',')]
        for i, val in enumerate(bs_vals):
            if val == 1: criteria_met.append(f"BS{i+1}")
    bp_match = re.search(r'BP=\[([^\]]+)\]', intervar_field)
    if bp_match:
        bp_vals = [int(x.strip()) for x in bp_match.group(1).split(',')]
        for i, val in enumerate(bp_vals):
            if val == 1: criteria_met.append(f"BP{i+1}")
    return (classification, criteria_met)

def parse_intervar_predictions(filepath: Path) -> tuple[dict[str, dict], dict[str, dict]]:
    predictions = {}
    gene_pos_index = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row.get('#Chr', '')
            if not chrom.startswith('chr'): chrom = f"chr{chrom}"
            pos = int(row.get('Start', 0))
            ref = row.get('Ref', '')
            alt = row.get('Alt', '')
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
            intervar_field = row.get(' InterVar: InterVar and Evidence ', '')
            classification, criteria_met = parse_intervar_classification(intervar_field)
            gene = row.get('Ref.Gene', '')
            pred_data = {
                "classification": classification, "criteria_met": criteria_met,
                "gene": gene, "consequence": row.get('Func.refGene', ''),
                "exonic_func": row.get('ExonicFunc.refGene', ''), "variant_key": variant_key
            }
            if variant_key not in predictions: predictions[variant_key] = []
            predictions[variant_key].append(pred_data)
            if gene:
                gene_pos_key = f"{gene.lower()}|{chrom}:{pos}"
                gene_pos_index[gene_pos_key] = pred_data
    return predictions, gene_pos_index

def parse_bias_predictions(filepath: Path) -> tuple[dict[str, dict], dict[str, dict], dict[str, dict], dict[str, dict], dict[str, dict]]:
    predictions = {}
    hgvsc_index = {}
    gene_hgvsc_index = {}
    hgvsg_index = {}
    position_index = defaultdict(list)
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row.get('chromosome', '')
            # ... (rest of logic) ...
            if not chrom.startswith('chr'): chrom = f"chr{chrom}"
            pos = int(row.get('position', 0))
            ref = row.get('refAllele', '')
            alt = row.get('altAllele', '')
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
            classification = row.get('acmgClassification', '').strip()
            classification_lower = classification.lower()
            if "pathogenic" in classification_lower and "likely" not in classification_lower: classification = "Pathogenic"
            elif "likely_pathogenic" in classification_lower or "likely pathogenic" in classification_lower: classification = "Likely Pathogenic"
            elif "uncertain" in classification_lower: classification = "Uncertain Significance"
            elif "likely_benign" in classification_lower or "likely benign" in classification_lower: classification = "Likely Benign"
            elif "benign" in classification_lower: classification = "Benign"
            else: classification = classification.title()
            criteria_str = row.get('acmgCriteria', '')
            criteria_met = [c.strip() for c in criteria_str.split(',') if c.strip()] if criteria_str else []
            gene = row.get('geneName', '')
            hgvsc = row.get('hgvsc', '')
            hgvsg = row.get('hgvsg', '')
            pred_data = {
                "classification": classification, "criteria_met": criteria_met,
                "gene": gene, "consequence": row.get('consequence', ''),
                "hgvsg": hgvsg, "hgvsc": hgvsc, "variant_key": variant_key,
                "ref": ref, "alt": alt, "chromosome": chrom, "position": pos # Add explicit fields for easier access
            }
            predictions[variant_key] = pred_data
            if hgvsc:
                hgvsc_index[hgvsc] = pred_data
                c_match = re.search(r'(c\.[^\s]+)', hgvsc)
                if c_match and gene:
                    gene_hgvsc_key = f"{gene.lower()}|{c_match.group(1)}"
                    gene_hgvsc_index[gene_hgvsc_key] = pred_data
            if hgvsg: hgvsg_index[hgvsg] = pred_data
            pos_key = f"{chrom}:{pos}"
            position_index[pos_key].append(pred_data)
    return predictions, hgvsc_index, gene_hgvsc_index, hgvsg_index, position_index


def normalize_indel_reproduce(ref: str, alt: str) -> tuple[str, str, int]:
    if not ref or not alt: return ref, alt, 0
    if ref[0] == alt[0]:
        new_ref = ref[1:] if len(ref) > 1 else '-'
        new_alt = alt[1:] if len(alt) > 1 else '-'
        offset = 1
        while new_ref and new_alt and new_ref != '-' and new_alt != '-':
            if len(new_ref) > 0 and len(new_alt) > 0 and new_ref[0] == new_alt[0]:
                new_ref = new_ref[1:] if len(new_ref) > 1 else '-'
                new_alt = new_alt[1:] if len(new_alt) > 1 else '-'
                offset += 1
            else:
                break
        return new_ref, new_alt, offset
    return ref, alt, 0

def parse_vcf_reproduce(vcf_path: Path) -> tuple[dict, dict]:
    variant_to_erepo = {}
    erepo_to_vcf = {}
    if not vcf_path.exists(): return {}, {}
    with open(vcf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            chrom, pos_str, vcf_id, ref, alt = parts[:5]
            pos = int(pos_str)
            info = parts[7]
            info_dict = dict((item.split('=') for item in info.split(';') if '=' in item))
            erepo_line = int(info_dict.get('EreppoLine', 0))
            if erepo_line == 0: continue
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
            variant_to_erepo[variant_key] = erepo_line
            erepo_to_vcf[erepo_line] = {
                'variant_key': variant_key, 'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                'vcf_id': vcf_id, 'variant_type': 'SNV'
            }
            new_ref, new_alt, offset = normalize_indel_reproduce(ref, alt)
            if offset > 0:
                normalized_key = f"{chrom}:{pos + offset}:{new_ref}:{new_alt}"
                variant_to_erepo[normalized_key] = erepo_line
    return variant_to_erepo, erepo_to_vcf
