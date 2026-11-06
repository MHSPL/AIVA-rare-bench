"""
Part 1 of the benchmark analysis: gene and phenotype extraction, and disease
categorization based on ClinGen Expert Panels. Writes the per-variant, per-gene
and per-category TSV summaries.
"""

import csv
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.categories import CATEGORY_ORDER, get_category

CLASS_COLUMNS = ["pathogenic", "likely_pathogenic", "vus", "likely_benign", "benign"]
CLASS_LABELS = ["Pathogenic", "Likely Pathogenic", "Uncertain Significance", "Likely Benign", "Benign"]


def load_jsonl(filepath: Path) -> list[dict]:
    """Load JSONL file with error handling."""
    results = []
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    try:
                        results.append(json.loads(line))
                    except json.JSONDecodeError as e:
                        print(f"Warning: Failed to parse JSON on line {line_num} in {filepath}: {e}")
                        continue
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}")
        raise
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        raise
    return results


def extract_gene_phenotype_data(filepaths: list[Path]) -> dict:
    """Extract unique gene-phenotype-variant combinations."""
    variants = {}
    for fp in filepaths:
        for r in load_jsonl(fp):
            vid = r.get("variant_id")
            if not vid or vid in variants:
                continue
            variants[vid] = {
                "variant_id": vid,
                "variant": r.get("variant") or "",
                "gene": r.get("gene", "Unknown"),
                "disease": r.get("disease", ""),
                "expert_panel": r.get("expert_panel", ""),
                "ground_truth": r.get("ground_truth", ""),
                "ground_truth_criteria": ", ".join(r.get("ground_truth_criteria", []))
            }
    return variants


def generate_gene_summary(variants: dict) -> dict:
    """Generate summary statistics by gene."""
    gene_data = defaultdict(lambda: {
        "count": 0,
        "phenotypes": set(),
        "expert_panels": set(),
        "classifications": defaultdict(int)
    })

    for v in variants.values():
        gene = v["gene"] or "Unknown"
        gene_data[gene]["count"] += 1
        if v["disease"]:
            gene_data[gene]["phenotypes"].add(v["disease"])
        if v["expert_panel"]:
            gene_data[gene]["expert_panels"].add(v["expert_panel"])
        if v["ground_truth"]:
            gene_data[gene]["classifications"][v["ground_truth"]] += 1

    for gene in gene_data:
        gene_data[gene]["phenotypes"] = sorted(gene_data[gene]["phenotypes"])
        gene_data[gene]["expert_panels"] = sorted(gene_data[gene]["expert_panels"])
        gene_data[gene]["classifications"] = dict(gene_data[gene]["classifications"])

    return dict(gene_data)


def write_gene_phenotype_files(variants: dict, gene_data: dict, output_dir: Path):
    """Write gene_phenotype_variants.tsv and gene_summary.tsv."""
    # Per-variant file
    variant_path = output_dir / "gene_phenotype_variants.tsv"
    fieldnames = ["variant_id", "variant", "gene", "disease", "expert_panel", "ground_truth", "ground_truth_criteria"]

    with open(variant_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(variants.values(), key=lambda x: (x["gene"] or "ZZZ", x["variant"])))

    print(f"  [1/8] Wrote {len(variants)} variants to {variant_path.name}")

    # Gene summary file
    gene_path = output_dir / "gene_summary.tsv"
    fieldnames = ["gene", "variant_count", "phenotypes", "expert_panels"] + CLASS_COLUMNS

    with open(gene_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for gene, data in sorted(gene_data.items(), key=lambda x: -x[1]["count"]):
            row = {
                "gene": gene,
                "variant_count": data["count"],
                "phenotypes": "; ".join(data["phenotypes"]),
                "expert_panels": "; ".join(data["expert_panels"]),
            }
            for col, label in zip(CLASS_COLUMNS, CLASS_LABELS):
                row[col] = data["classifications"].get(label, 0)
            writer.writerow(row)

    print(f"  [2/8] Wrote {len(gene_data)} genes to {gene_path.name}")


def categorize_and_write(variants: dict, output_dir: Path):
    """Categorize variants and write phenotype_categories.tsv and category_summary.tsv."""
    for v in variants.values():
        v["disease_category"] = get_category(v["expert_panel"])

    # Write phenotype_categories.tsv
    cat_path = output_dir / "phenotype_categories.tsv"
    fieldnames = list(next(iter(variants.values())).keys())

    with open(cat_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants.values())

    print(f"  [3/8] Wrote {len(variants)} categorized variants to {cat_path.name}")

    # Category summary
    category_counts = Counter(v["disease_category"] for v in variants.values())
    gene_by_category = defaultdict(set)
    vcep_by_category = defaultdict(set)
    classification_by_category = defaultdict(lambda: defaultdict(int))

    for v in variants.values():
        cat = v["disease_category"]
        gene_by_category[cat].add(v["gene"])
        vcep_by_category[cat].add(v["expert_panel"])
        classification_by_category[cat][v["ground_truth"]] += 1

    summary_path = output_dir / "category_summary.tsv"

    def sort_key(cat):
        try:
            return (CATEGORY_ORDER.index(cat), -category_counts[cat])
        except ValueError:
            return (len(CATEGORY_ORDER), -category_counts[cat])

    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "variant_count", "vcep_count", "gene_count"]
                        + CLASS_COLUMNS + ["expert_panels", "genes"])

        for cat in sorted(category_counts.keys(), key=sort_key):
            cls_counts = classification_by_category[cat]
            writer.writerow(
                [cat, category_counts[cat], len(vcep_by_category[cat]), len(gene_by_category[cat])]
                + [cls_counts.get(label, 0) for label in CLASS_LABELS]
                + ["; ".join(sorted(v for v in vcep_by_category[cat] if v)),
                   "; ".join(sorted(gene_by_category[cat]))]
            )

    print(f"  [4/8] Wrote category summary to {summary_path.name}")
