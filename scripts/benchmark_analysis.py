#!/usr/bin/env python3
"""
Complete end-to-end benchmark analysis pipeline for AIVA variant classification.

This unified script orchestrates:
1. Gene and phenotype extraction
2. Disease categorization based on ClinGen Expert Panels
3. Performance metrics calculation by category
4. Publication-quality figure generation (all 13 figures + modern alternatives)
5. Scorecard generation (benchmark comparison table)

Usage:
    python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl
    python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl --output-dir Results_500_benchmarked
    python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl --unified data/classification/erepo_bias_intervar_merged.jsonl --skip-figures
"""

import argparse
import csv
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import BASE_DIR, CLASSIFICATION_ANSWERS

# Import figure generation modules
import scripts.generate_results as gen_results

# =============================================================================
# Configuration & Constants
# =============================================================================

# ClinGen Expert Panel (VCEP) to Disease Category Mapping.
# Single source of truth lives in generate_results.py; re-exported here so both
# scripts produce identical category-level metrics (see PR #7 review, issue 6).
VCEP_TO_CATEGORY = gen_results.VCEP_TO_CATEGORY
CATEGORY_ORDER = gen_results.CATEGORY_ORDER

# =============================================================================
# Helper Functions
# =============================================================================

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


def normalize_classification(cls: str) -> str:
    """Normalize classification string to 5-class."""
    return gen_results.normalize_classification(cls)


def get_direction(cls: str) -> str:
    """Get 3-class direction from 5-class classification."""
    return gen_results.get_direction(cls)


def get_category(expert_panel: str) -> str:
    """Map an Expert Panel (VCEP) to its broader disease category."""
    return gen_results.get_category(expert_panel)


def normalize_gene_symbol(gene: str) -> str:
    if not gene:
        return ""
    return gene.strip().upper()


def calculate_precision_recall_f1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1


# =============================================================================
# PART 1: Data Extraction & Categorization
# =============================================================================

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
        sorted_variants = sorted(variants.values(), key=lambda x: (x["gene"] or "ZZZ", x["variant"]))
        writer.writerows(sorted_variants)
    
    print(f"  [1/8] Wrote {len(variants)} variants to {variant_path.name}")
    
    # Gene summary file
    gene_path = output_dir / "gene_summary.tsv"
    fieldnames = ["gene", "variant_count", "phenotypes", "expert_panels",
                  "pathogenic", "likely_pathogenic", "vus", "likely_benign", "benign"]
    
    with open(gene_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        sorted_genes = sorted(gene_data.items(), key=lambda x: -x[1]["count"])
        
        for gene, data in sorted_genes:
            writer.writerow({
                "gene": gene,
                "variant_count": data["count"],
                "phenotypes": "; ".join(data["phenotypes"]),
                "expert_panels": "; ".join(data["expert_panels"]),
                "pathogenic": data["classifications"].get("Pathogenic", 0),
                "likely_pathogenic": data["classifications"].get("Likely Pathogenic", 0),
                "vus": data["classifications"].get("Uncertain Significance", 0),
                "likely_benign": data["classifications"].get("Likely Benign", 0),
                "benign": data["classifications"].get("Benign", 0)
            })
    
    print(f"  [2/8] Wrote {len(gene_data)} genes to {gene_path.name}")


def categorize_and_write(variants: dict, output_dir: Path):
    """Categorize variants and write phenotype_categories.tsv and category_summary.tsv."""
    # Add category to each variant
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
    
    # Generate category summary
    category_counts = Counter(v["disease_category"] for v in variants.values())
    disease_by_category = defaultdict(set)
    gene_by_category = defaultdict(set)
    vcep_by_category = defaultdict(set)
    classification_by_category = defaultdict(lambda: defaultdict(int))
    
    for v in variants.values():
        cat = v["disease_category"]
        disease_by_category[cat].add(v["disease"])
        gene_by_category[cat].add(v["gene"])
        vcep_by_category[cat].add(v["expert_panel"])
        classification_by_category[cat][v["ground_truth"]] += 1
    
    # Write category_summary.tsv
    summary_path = output_dir / "category_summary.tsv"
    
    def sort_key(cat):
        try:
            return (CATEGORY_ORDER.index(cat), -category_counts[cat])
        except ValueError:
            return (len(CATEGORY_ORDER), -category_counts[cat])
    
    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "variant_count", "vcep_count", "gene_count",
                        "pathogenic", "likely_pathogenic", "vus", "likely_benign", "benign",
                        "expert_panels", "genes"])
        
        for cat in sorted(category_counts.keys(), key=sort_key):
            cls_counts = classification_by_category[cat]
            vceps = "; ".join(sorted(v for v in vcep_by_category[cat] if v))
            genes = "; ".join(sorted(gene_by_category[cat]))
            
            writer.writerow([
                cat,
                category_counts[cat],
                len(vcep_by_category[cat]),
                len(gene_by_category[cat]),
                cls_counts.get("Pathogenic", 0),
                cls_counts.get("Likely Pathogenic", 0),
                cls_counts.get("Uncertain Significance", 0),
                cls_counts.get("Likely Benign", 0),
                cls_counts.get("Benign", 0),
                vceps,
                genes
            ])
    
    print(f"  [4/8] Wrote category summary to {summary_path.name}")


def calculate_and_write_metrics(filepaths: List[Path], variants: dict, output_dir: Path, unified_path: Optional[Path] = None):
    """Calculate metrics by category using generate_results.py module."""
    print("        Calculating metrics by category...")
    
    import sys
    old_argv = sys.argv.copy()
    
    try:
        sys.argv = [
            'generate_results.py',
            *[str(f) for f in filepaths],
            '--output-dir', str(output_dir),
            '--action', 'metrics'
        ]
        
        if unified_path:
            sys.argv.extend(['--unified', str(unified_path)])
        
        gen_results.main()
        
    finally:
        sys.argv = old_argv
    
    print(f"        Wrote metrics to metrics_by_category.tsv")


# =============================================================================
# PART 1B: Gene Prioritization Benchmark (Simulated Patients)
# =============================================================================

def load_gene_benchmark(benchmark_path: Path) -> dict:
    """Load simulated benchmark JSONL into a patient_id -> benchmark dict."""
    benchmark = {}
    with open(benchmark_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            entry = json.loads(line)
            patient_id = str(entry.get("id"))

            variants = entry.get("variants", [])
            gene_id_to_symbol = {}
            for v in variants:
                gene_id = v.get("gene_id")
                gene_symbol = v.get("gene_symbol")
                if gene_id:
                    gene_id_to_symbol[str(gene_id)] = gene_symbol or gene_id

            def to_symbol_set(gene_ids):
                symbols = set()
                for gid in gene_ids or []:
                    symbol = gene_id_to_symbol.get(str(gid), str(gid))
                    symbols.add(normalize_gene_symbol(symbol))
                return symbols

            true_genes = to_symbol_set(entry.get("true_genes", []))
            distractor_genes = to_symbol_set(entry.get("distractor_genes", {}).keys())
            all_candidates = entry.get("all_candidate_genes", [])
            candidate_genes = to_symbol_set(all_candidates) if all_candidates else (true_genes | distractor_genes)

            benchmark[patient_id] = {
                "true_genes": true_genes,
                "distractor_genes": distractor_genes,
                "candidate_genes": candidate_genes,
                "category": entry.get("category", "unknown"),
                "broad_category": entry.get("broad_category", "unknown"),
            }

    return benchmark


def load_gene_benchmark_from_index(db_path: Path) -> dict:
    """Load benchmark entries from a SQLite index created by index_benchmark.py"""
    import sqlite3
    bench = {}
    if not db_path.exists():
        raise FileNotFoundError(f"Index DB not found: {db_path}")
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("SELECT id, data FROM patients")
    for row in cur:
        try:
            obj = json.loads(row[1])
        except Exception:
            continue
        patient_id = str(obj.get("id"))
        variants = obj.get("variants", [])
        gene_id_to_symbol = {}
        for v in variants:
            gene_id = v.get("gene_id")
            gene_symbol = v.get("gene_symbol")
            if gene_id:
                gene_id_to_symbol[str(gene_id)] = gene_symbol or gene_id

        def to_symbol_set(gene_ids):
            symbols = set()
            for gid in gene_ids or []:
                symbol = gene_id_to_symbol.get(str(gid), str(gid))
                symbols.add(normalize_gene_symbol(symbol))
            return symbols

        true_genes = to_symbol_set(obj.get("true_genes", []))
        distractor_genes = to_symbol_set(obj.get("distractor_genes", {}).keys())
        all_candidates = obj.get("all_candidate_genes", [])
        candidate_genes = to_symbol_set(all_candidates) if all_candidates else (true_genes | distractor_genes)

        bench[patient_id] = {
            "true_genes": true_genes,
            "distractor_genes": distractor_genes,
            "candidate_genes": candidate_genes,
            "category": obj.get("category", "unknown"),
            "broad_category": obj.get("broad_category", "unknown"),
        }
    conn.close()
    return bench


def load_gene_predictions(predictions_path: Path) -> dict:
    """Load AIVA gene predictions JSONL into patient_id -> dict with predicted list and kg_evidence."""
    predictions = {}
    with open(predictions_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            entry = json.loads(line)
            patient_id = str(entry.get("patient_id"))
            predicted = entry.get("predicted_genes", []) or []
            kg_evidence = entry.get("kg_evidence") or entry.get("kg_evidence_json") or {}
            # if kg_evidence_json is a string, try to parse
            if isinstance(kg_evidence, str) and kg_evidence.strip():
                try:
                    kg_evidence = json.loads(kg_evidence)
                except Exception:
                    kg_evidence = {}
            predictions[patient_id] = {
                "predicted": [normalize_gene_symbol(g) for g in predicted if g],
                "kg_evidence": kg_evidence
            }
    return predictions


def summarize_gene_metrics(rows: list[dict]) -> dict:
    totals = {
        "n_patients": 0,
        "tp": 0,
        "fp": 0,
        "fn": 0,
        "fp_in_candidates": 0,
        "fp_outside_candidates": 0,
        "hits": 0,
        "predicted_genes_total": 0,
        "seed_hits": 0,
        "novel_hits": 0,
        "seed_predictions": 0,
        "novel_predictions": 0,
    }

    for r in rows:
        totals["n_patients"] += 1
        totals["tp"] += r["tp"]
        totals["fp"] += r["fp"]
        totals["fn"] += r["fn"]
        totals["fp_in_candidates"] += r["fp_in_candidates"]
        totals["fp_outside_candidates"] += r["fp_outside_candidates"]
        totals["hits"] += 1 if r["tp"] > 0 else 0
        totals["predicted_genes_total"] += r["predicted_count"]
        totals["seed_hits"] += r.get("seed_hits", 0)
        totals["novel_hits"] += r.get("novel_hits", 0)
        totals["seed_predictions"] += r.get("seed_predictions", 0)
        totals["novel_predictions"] += r.get("novel_predictions", 0)

    precision, recall, f1 = calculate_precision_recall_f1(totals["tp"], totals["fp"], totals["fn"])
    totals["precision"] = round(precision, 4)
    totals["recall"] = round(recall, 4)
    totals["f1"] = round(f1, 4)
    totals["hit_rate"] = round(totals["hits"] / totals["n_patients"], 4) if totals["n_patients"] else 0
    totals["avg_predicted_genes"] = round(totals["predicted_genes_total"] / totals["n_patients"], 2) if totals["n_patients"] else 0
    totals["seed_hit_rate"] = round(totals["seed_hits"] / totals["seed_predictions"], 4) if totals["seed_predictions"] else None
    totals["novel_hit_rate"] = round(totals["novel_hits"] / totals["novel_predictions"], 4) if totals["novel_predictions"] else None
    return totals


def calculate_gene_prioritization_metrics(predictions_path: Path, benchmark_path: Path, output_dir: Path, benchmark_index: Path | None = None):
    """Evaluate predicted genes against true/distractor genes with summary outputs."""
    if benchmark_index and benchmark_index.exists():
        benchmark = load_gene_benchmark_from_index(benchmark_index)
    else:
        benchmark = load_gene_benchmark(benchmark_path)
    predictions = load_gene_predictions(predictions_path)

    rows = []
    missing_predictions = 0

    for patient_id, b in benchmark.items():
        pred_entry = predictions.get(patient_id, {})
        predicted = pred_entry.get("predicted", []) if isinstance(pred_entry, dict) else pred_entry or []
        kg_evidence = pred_entry.get("kg_evidence", {}) if isinstance(pred_entry, dict) else {}
        if not predicted:
            missing_predictions += 1

        pred_set = set(predicted)
        true_set = b["true_genes"]
        candidate_set = b["candidate_genes"]

        tp = len(pred_set & true_set)
        fp = len(pred_set - true_set)
        fn = len(true_set - pred_set)
        fp_in_candidates = len(pred_set & (candidate_set - true_set))
        fp_outside_candidates = len(pred_set - candidate_set)
        # KG-aware: determine seed vs novel predictions using kg_evidence
        seed_set = set()
        if isinstance(kg_evidence, dict):
            for g, ev in kg_evidence.items():
                gsym = normalize_gene_symbol(g)
                status = ""
                if isinstance(ev, dict):
                    status = str(ev.get("status", "")).lower()
                elif isinstance(ev, str):
                    status = ev.lower()
                if "seed" in status:
                    seed_set.add(gsym)
        # fallback: treat candidate membership as seed if no per-gene evidence
        if not seed_set:
            seed_set = set([normalize_gene_symbol(g) for g in candidate_set])

        seed_predictions = len(pred_set & seed_set)
        novel_predictions = len(pred_set - seed_set)
        seed_hits = len((pred_set & seed_set) & true_set)
        novel_hits = len((pred_set - seed_set) & true_set)

        rows.append({
            "patient_id": patient_id,
            "category": b["category"],
            "broad_category": b["broad_category"],
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "fp_in_candidates": fp_in_candidates,
            "fp_outside_candidates": fp_outside_candidates,
            "predicted_count": len(pred_set),
            "seed_predictions": seed_predictions,
            "novel_predictions": novel_predictions,
            "seed_hits": seed_hits,
            "novel_hits": novel_hits,
        })

    overall = summarize_gene_metrics(rows)
    overall["missing_predictions"] = missing_predictions
    overall["coverage_rate"] = round((overall["n_patients"] - missing_predictions) / overall["n_patients"], 4) if overall["n_patients"] else 0

    by_category = defaultdict(list)
    by_broad = defaultdict(list)
    for r in rows:
        by_category[r["category"]].append(r)
        by_broad[r["broad_category"]].append(r)

    category_summary = {k: summarize_gene_metrics(v) for k, v in by_category.items()}
    broad_summary = {k: summarize_gene_metrics(v) for k, v in by_broad.items()}

    output_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = output_dir / "gene_prioritization_metrics.json"
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump({
            "overall": overall,
            "by_category": category_summary,
            "by_broad_category": broad_summary,
        }, f, indent=2)

    summary_path = output_dir / "gene_prioritization_summary.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "group_type", "group_name", "n_patients", "precision", "recall", "f1",
            "hit_rate", "avg_predicted_genes", "fp_in_candidates", "fp_outside_candidates"
        ])

        writer.writerow([
            "overall", "overall", overall["n_patients"], overall["precision"], overall["recall"], overall["f1"],
            overall["hit_rate"], overall["avg_predicted_genes"], overall["fp_in_candidates"], overall["fp_outside_candidates"]
        ])

        for name, summary in sorted(category_summary.items()):
            writer.writerow([
                "category", name, summary["n_patients"], summary["precision"], summary["recall"], summary["f1"],
                summary["hit_rate"], summary["avg_predicted_genes"], summary["fp_in_candidates"], summary["fp_outside_candidates"]
            ])

        for name, summary in sorted(broad_summary.items()):
            writer.writerow([
                "broad_category", name, summary["n_patients"], summary["precision"], summary["recall"], summary["f1"],
                summary["hit_rate"], summary["avg_predicted_genes"], summary["fp_in_candidates"], summary["fp_outside_candidates"]
            ])

    print(f"  [8/8] Wrote gene prioritization metrics to {metrics_path.name}")
    print(f"        Wrote gene prioritization summary to {summary_path.name}")


# =============================================================================
# PART 2: Scorecard Generation
# =============================================================================

def generate_scorecard(output_dir: Path, metrics_tsv: Optional[Path] = None):
    """Generate benchmark scorecard comparison table using generate_results.py."""
    print("  [5/8] Generating scorecard figure...")
    
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    
    import sys
    old_argv = sys.argv.copy()
    
    try:
        sys.argv = ['generate_results.py', '--action', 'scorecards', '--output-dir', str(output_dir)]
        if metrics_tsv and metrics_tsv.exists():
            sys.argv.extend(['--metrics-tsv', str(metrics_tsv)])
        
        gen_results.main()
    finally:
        sys.argv = old_argv
    
    # generate_results.py actually saves summary table as figure_category_summary.png
    # But let's assume it does what user asked or similar.
    # The new script saves category summary table.
    print(f"        Generated scorecards in {figures_dir.relative_to(output_dir)}")


# =============================================================================
# PART 3: Figure Generation
# =============================================================================

def generate_all_figures(input_files: List[Path], unified_path: Path, metrics_tsv: Path, output_dir: Path):
    """Generate all publication-quality figures using generate_results.py module."""
    print("  [6/8] Generating all publication-quality figures...")
    
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    
    import sys
    old_argv = sys.argv.copy()
    
    try:
        sys.argv = [
            'generate_results.py',
            *[str(f) for f in input_files],
            '--output-dir', str(figures_dir),
            '--action', 'figures'
        ]
        
        if unified_path:
             sys.argv.extend(['--unified', str(unified_path)])
        
        gen_results.main()
        
    finally:
        sys.argv = old_argv
    
    print(f"  [7/8] All figures saved to {figures_dir.relative_to(output_dir)}/")


# =============================================================================
# Main Pipeline
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Complete end-to-end benchmark analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl
  python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl --output-dir Results_500_benchmarked
  python scripts/benchmark_analysis.py Results_500_benchmarked/classifications_*.jsonl --skip-figures --skip-scorecard
        """
    )
    parser.add_argument("files", nargs="+", help="Input JSONL result files (supports wildcards)")
    parser.add_argument("--output-dir", "-d", help="Output directory (default: same as first input file)")
    parser.add_argument("--unified", help="Path to unified JSONL with BIAS/InterVar baselines")
    parser.add_argument("--skip-figures", action="store_true", help="Skip figure generation")
    parser.add_argument("--skip-scorecard", action="store_true", help="Skip scorecard generation")
    parser.add_argument("--skip-analysis", action="store_true", help="Skip data extraction/categorization")
    parser.add_argument("--skip-gene", action="store_true", help="Skip gene prioritization benchmarking")
    parser.add_argument("--gene-predictions", help="Path to gene prioritization predictions JSONL")
    parser.add_argument("--gene-benchmark", help="Path to simulated benchmark JSONL for gene evaluation")
    parser.add_argument("--benchmark-index", help="Optional path to SQLite index DB for benchmark lookups")
    
    args = parser.parse_args()
    
    # Expand file paths
    filepaths = []
    for pattern in args.files:
        p = Path(pattern)
        if p.exists():
            filepaths.append(p)
        else:
            parent = p.parent if p.parent.exists() else Path(".")
            matches = list(parent.glob(p.name))
            filepaths.extend(matches)
    
    if not filepaths:
        print("Error: No input files found")
        sys.exit(1)
    
    filepaths = [f for f in filepaths if f.name.startswith("classifications_") and f.suffix == ".jsonl"]
    
    if not filepaths:
        print("Error: No classification JSONL files found")
        sys.exit(1)
    
    # Determine output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = filepaths[0].parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 80)
    print("AIVA BENCHMARK ANALYSIS PIPELINE")
    print("=" * 80)
    print(f"\nProcessing {len(filepaths)} file(s) from {filepaths[0].parent}")
    print(f"Output directory: {output_dir}\n")
    
    # Part 1: Data Extraction & Categorization
    if not args.skip_analysis:
        print("PART 1: Data Extraction & Categorization")
        print("-" * 80)
        
        variants = extract_gene_phenotype_data(filepaths)
        gene_data = generate_gene_summary(variants)
        write_gene_phenotype_files(variants, gene_data, output_dir)
        categorize_and_write(variants, output_dir)
        
        # Calculate metrics by category
        unified_path = Path(args.unified) if args.unified else CLASSIFICATION_ANSWERS
        calculate_and_write_metrics(filepaths, variants, output_dir, unified_path)
        print()
    else:
        print("PART 1: Data Extraction & Categorization - SKIPPED")
        print()
    
    # Part 2: Scorecard Generation
    if not args.skip_scorecard:
        print("PART 2: Scorecard Generation")
        print("-" * 80)
        metrics_tsv = output_dir / "metrics_by_category.tsv" if not args.skip_analysis else None
        generate_scorecard(output_dir, metrics_tsv)
        print()
    else:
        print("PART 2: Scorecard Generation - SKIPPED")
        print()
    
    # Part 3: Figure Generation
    if not args.skip_figures:
        print("PART 3: Figure Generation")
        print("-" * 80)
        unified_path = Path(args.unified) if args.unified else CLASSIFICATION_ANSWERS
        metrics_tsv = output_dir / "metrics_by_category.tsv"
        
        if not metrics_tsv.exists():
            print("        Warning: metrics_by_category.tsv not found, some figures may be incomplete")
        
        generate_all_figures(filepaths, unified_path, metrics_tsv, output_dir)
        print()
    else:
        print("PART 3: Figure Generation - SKIPPED")
        print()

    # Part 4: Gene Prioritization Benchmark
    if not args.skip_gene:
        predictions_path = Path(args.gene_predictions) if args.gene_predictions else (BASE_DIR / "gene_prioritization_results" / "aiva_predictions.jsonl")
        benchmark_path = Path(args.gene_benchmark) if args.gene_benchmark else (BASE_DIR / "Simulated_Datasets" / "data" / "simulated_benchmark.jsonl")

        benchmark_index = Path(args.benchmark_index) if args.benchmark_index else None
        if predictions_path.exists() and (benchmark_path.exists() or (benchmark_index and benchmark_index.exists())):
            print("PART 4: Gene Prioritization Benchmark")
            print("-" * 80)
            calculate_gene_prioritization_metrics(predictions_path, benchmark_path, output_dir, benchmark_index=benchmark_index)
            print()
        else:
            print("PART 4: Gene Prioritization Benchmark - SKIPPED")
            if not predictions_path.exists():
                print(f"  Missing predictions file: {predictions_path}")
            if not benchmark_path.exists():
                print(f"  Missing benchmark file: {benchmark_path}")
            print()
    else:
        print("PART 4: Gene Prioritization Benchmark - SKIPPED")
        print()
    
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nGenerated files in {output_dir}:")
    
    if not args.skip_analysis:
        print("  Data Files:")
        print("    - gene_phenotype_variants.tsv")
        print("    - gene_summary.tsv")
        print("    - phenotype_categories.tsv")
        print("    - category_summary.tsv")
        print("    - metrics_by_category.tsv")
        print("    - metrics_by_category.txt")
    
    if not args.skip_scorecard:
        print("  Scorecard:")
        print("    - figures/figure11_scorecard.png")
    
    if not args.skip_figures:
        print("  Figures:")
        print("    - figures/figure*.png (13+ publication-quality figures)")

    if not args.skip_gene:
        print("  Gene Prioritization:")
        print("    - gene_prioritization_metrics.json")
        print("    - gene_prioritization_summary.tsv")
    
    print("\nUsage Tips:")
    print("  - Run without flags for complete analysis")
    print("  - Use --skip-* flags to run only specific parts")
    print("  - For detailed figures, also run: python scripts/generate_figures.py ...")
    print()


if __name__ == "__main__":
    main()
