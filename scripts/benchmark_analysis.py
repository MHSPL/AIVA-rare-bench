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
import sys
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import BASE_DIR, CLASSIFICATION_ANSWERS

import scripts.generate_results as gen_results
from scripts.categories import VCEP_TO_CATEGORY, CATEGORY_ORDER, get_category  # noqa: F401  (re-exported)
from scripts.metrics import (  # noqa: F401  (re-exported)
    normalize_classification, get_direction, calculate_precision_recall_f1,
)
from scripts.variant_summaries import (  # noqa: F401  (re-exported)
    load_jsonl, extract_gene_phenotype_data, generate_gene_summary,
    write_gene_phenotype_files, categorize_and_write,
)
from scripts.gene_prioritization import (  # noqa: F401  (re-exported)
    normalize_gene_symbol, load_gene_benchmark, load_gene_benchmark_from_index,
    load_gene_predictions, summarize_gene_metrics, calculate_gene_prioritization_metrics,
)


def run_generate_results(argv: list[str]):
    """Run generate_results.main() with the given command line arguments."""
    old_argv = sys.argv.copy()
    try:
        sys.argv = ['generate_results.py', *argv]
        gen_results.main()
    finally:
        sys.argv = old_argv


# =============================================================================
# PART 1: Data Extraction & Categorization (see scripts/variant_summaries.py)
# =============================================================================

def calculate_and_write_metrics(filepaths: List[Path], variants: dict, output_dir: Path, unified_path: Optional[Path] = None):
    """Calculate metrics by category using generate_results.py module."""
    print("        Calculating metrics by category...")
    argv = [*[str(f) for f in filepaths], '--output-dir', str(output_dir), '--action', 'metrics']
    if unified_path:
        argv.extend(['--unified', str(unified_path)])
    run_generate_results(argv)
    print(f"        Wrote metrics to metrics_by_category.tsv")


# =============================================================================
# PART 2: Scorecard Generation
# =============================================================================

def generate_scorecard(output_dir: Path, metrics_tsv: Optional[Path] = None):
    """Generate benchmark scorecard comparison table using generate_results.py."""
    print("  [5/8] Generating scorecard figure...")
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    argv = ['--action', 'scorecards', '--output-dir', str(output_dir)]
    if metrics_tsv and metrics_tsv.exists():
        argv.extend(['--metrics-tsv', str(metrics_tsv)])
    run_generate_results(argv)
    print(f"        Generated scorecards in {figures_dir.relative_to(output_dir)}")


# =============================================================================
# PART 3: Figure Generation
# =============================================================================

def generate_all_figures(input_files: List[Path], unified_path: Path, metrics_tsv: Path, output_dir: Path):
    """Generate all publication-quality figures using generate_results.py module."""
    print("  [6/8] Generating all publication-quality figures...")
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    argv = [*[str(f) for f in input_files], '--output-dir', str(figures_dir), '--action', 'figures']
    if unified_path:
        argv.extend(['--unified', str(unified_path)])
    run_generate_results(argv)
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
