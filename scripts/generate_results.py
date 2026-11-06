#!/usr/bin/env python3
"""
Results generation for RareDx Bench.

Entry point that ties together:
1. Category-based metric calculation (scripts/category_metrics.py)
2. Figure generation: confusion matrices, error analysis (scripts/figures.py)
3. Scorecard tables (scripts/summary_tables.py)

Core metric utilities live in scripts/metrics.py and disease categories in
scripts/categories.py; both are re-exported here for existing importers.

Usage:
    # 1. Calculate metrics by category and save TSV
    python scripts/generate_results.py classifications.jsonl --action metrics --output-dir results

    # 2. Generate all figures and scorecards from results & TSV
    python scripts/generate_results.py classifications.jsonl --metrics-tsv results/metrics_by_category.tsv --action all

    # Importable as a library
    from scripts.generate_results import calculate_category_metrics, load_results
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import CLASSIFICATION_LABELS, CLASSIFICATION_TO_DIRECTION, DIRECTION_LABELS, CLASSIFICATION_ANSWERS
from scripts.categories import (  # noqa: F401  (re-exported)
    TOOL_SHORT_NAMES, VCEP_TO_CATEGORY, CATEGORY_ORDER, get_short_name, get_category,
)
from scripts.metrics import (  # noqa: F401  (re-exported)
    load_results, normalize_classification, get_direction, normalize_criterion, get_criterion_type,
    calculate_precision_recall_f1, calculate_class_metrics, calculate_direction_metrics,
    calculate_macro_f1, calculate_criteria_metrics, extract_ground_truth_and_predictions,
)
from scripts.category_metrics import (  # noqa: F401  (re-exported)
    calculate_specificity, calculate_full_metrics, calculate_category_metrics, save_metrics_tsv,
)

# Plotting needs numpy, pandas and matplotlib; the metrics action works without them.
try:
    import numpy as np
    import pandas as pd
    import matplotlib  # noqa: F401
except ImportError:
    np = None
    pd = None
else:
    from scripts.figures import plot_confusion_matrix, plot_error_analysis
    from scripts.summary_tables import (
        category_macro_f1_table, generate_category_summary_table, generate_combined_summary,
    )


def run_metrics(input_files, unified, out_dir: Path):
    """Score every model in the input files, plus the BIAS-2015 and InterVar baselines."""
    print("Calculating metrics by category...")
    all_variants = []
    models = set()
    variant_ids = set()
    for f in input_files:
        if not Path(f).exists():
            continue
        loaded = load_results(Path(f))
        all_variants.extend(loaded)
        for v in loaded:
            if v.get("model_used"):
                models.add(v["model_used"])
            if v.get("variant_id"):
                variant_ids.add(v["variant_id"])

    unified_path = Path(unified) if unified else CLASSIFICATION_ANSWERS
    all_metrics = {m: calculate_category_metrics(all_variants, m) for m in models}

    if unified_path.exists():
        unified_filtered = [v for v in load_results(unified_path) if v.get("id") in variant_ids]
        if unified_filtered:
            all_metrics["BIAS-2015"] = calculate_category_metrics(unified_filtered, None, "bias_classification")
            all_metrics["InterVar"] = calculate_category_metrics(unified_filtered, None, "intervar_classification")

    save_metrics_tsv(all_metrics, out_dir / "metrics_by_category.tsv")


def _confusion(truth, preds, labels):
    cm = np.zeros((len(labels), len(labels)), dtype=int)
    for g, p in zip(truth, preds):
        if g in labels and p in labels:
            cm[labels.index(g), labels.index(p)] += 1
    return cm


def run_figures(input_files, out_dir: Path):
    """Confusion matrices over all inputs, and the error analysis of the best model."""
    print("Generating figures...")
    all_results = []
    for f in input_files:
        all_results.extend(load_results(Path(f)))
    gt, pred, gt_dir, pred_dir = extract_ground_truth_and_predictions(all_results)

    plot_confusion_matrix(_confusion(gt, pred, CLASSIFICATION_LABELS), CLASSIFICATION_LABELS, "5-Class", out_dir / "confusion_matrix_5class.png")
    plot_confusion_matrix(_confusion(gt_dir, pred_dir, DIRECTION_LABELS), DIRECTION_LABELS, "3-Class", out_dir / "confusion_matrix_3class.png")

    print("Generating Figure 10: Error analysis...")
    aiva_data_dict = {}
    for f in input_files:
        recs = load_results(Path(f))
        if recs:
            model = recs[0].get("model_used") or Path(f).stem
            aiva_data_dict[model] = recs
    plot_error_analysis(aiva_data_dict, out_dir)


def run_scorecards(metrics_tsv, out_dir: Path):
    """Category summary table and combined summary from metrics_by_category.tsv."""
    print("Generating scorecards...")
    df = pd.read_csv(metrics_tsv, sep='\t')
    generate_category_summary_table(category_macro_f1_table(df), out_dir)
    generate_combined_summary(metrics_tsv, out_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", nargs='*', help="JSONL files")
    parser.add_argument("--unified", help="Data path for baselines")
    parser.add_argument("--metrics-tsv", help="Path to TSV")
    parser.add_argument("--output-dir", default="figures")
    parser.add_argument("--action", choices=['all', 'metrics', 'figures', 'scorecards'], default='all')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.action in ['all', 'metrics']:
        run_metrics(args.input_files, args.unified, out_dir)
    if args.action in ['all', 'figures']:
        run_figures(args.input_files, out_dir)
    if args.action in ['all', 'scorecards'] and args.metrics_tsv:
        run_scorecards(args.metrics_tsv, out_dir)


if __name__ == "__main__":
    main()
