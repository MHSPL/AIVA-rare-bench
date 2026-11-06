"""
Per-category classification metrics for RareDx Bench.

Groups variants by disease category and scores directional (P/LP, VUS, B/LB)
and exact 5-class calls in each, then writes metrics_by_category.tsv.
"""

import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import CLASSIFICATION_LABELS
from scripts.categories import CATEGORY_ORDER, get_category
from scripts.metrics import calculate_precision_recall_f1, get_direction, normalize_classification

DIRECTION_TARGETS = ["pathogenic", "vus", "benign"]


def calculate_specificity(tn, fp):
    return tn / (tn + fp) if (tn + fp) > 0 else 0


def calculate_full_metrics(tp, fp, fn, tn):
    p, r, f1 = calculate_precision_recall_f1(tp, fp, fn)
    spec = calculate_specificity(tn, fp)
    return {"sensitivity": round(r, 4), "specificity": round(spec, 4), "precision": round(p, 4),
            "f1": round(f1, 4), "tp": tp, "fp": fp, "fn": fn, "tn": tn}


def _confusion_counts(truth, preds, target, n):
    tp = sum(1 for g, p in zip(truth, preds) if g == target and p == target)
    fp = sum(1 for g, p in zip(truth, preds) if g != target and p == target)
    fn = sum(1 for g, p in zip(truth, preds) if g == target and p != target)
    return tp, fp, fn, n - tp - fp - fn


def calculate_category_metrics(variants: List[Dict], model: str, pred_field: str = "aiva_classification") -> Dict:
    by_cat = defaultdict(list)
    for v in variants:
        if model and pred_field == "aiva_classification" and v.get("model_used") != model:
            continue
        by_cat[get_category(v.get("expert_panel", ""))].append(v)

    results = {}
    for cat, cat_variants in by_cat.items():
        pd_l, gt_l, pe_l, gte_l = [], [], [], []
        for v in cat_variants:
            p = v.get(pred_field, "")
            g = v.get("ground_truth", "") or v.get("classification", "")
            if p and g and not v.get("error"):
                pd_l.append(get_direction(p))
                gt_l.append(get_direction(g))
                pe_l.append(normalize_classification(p))
                gte_l.append(normalize_classification(g))

        n = len(pd_l)
        if n == 0:
            continue

        m = {"n_variants": n}
        # Directional metrics
        gt_low = [g.lower() for g in gt_l]
        pd_low = [p.lower() for p in pd_l]
        for target in DIRECTION_TARGETS:
            m[target] = calculate_full_metrics(*_confusion_counts(gt_low, pd_low, target, n))

        # Accuracy
        m["directional_accuracy"] = round(sum(1 for g, p in zip(gt_l, pd_l) if g == p) / n * 100, 2)
        m["exact_accuracy"] = round(sum(1 for g, p in zip(gte_l, pe_l) if g == p) / n * 100, 2)

        # Exact class metrics
        for cls in CLASSIFICATION_LABELS:
            m[f"exact_{cls}"] = calculate_full_metrics(*_confusion_counts(gte_l, pe_l, cls, n))

        results[cat] = m
    return results


TSV_HEADER = [
    "category", "model", "n_variants", "directional_accuracy",
    "plp_sensitivity", "plp_specificity", "plp_precision", "plp_f1",
    "vus_sensitivity", "vus_specificity", "vus_precision", "vus_f1",
    "blb_sensitivity", "blb_specificity", "blb_precision", "blb_f1",
    "exact_accuracy", "exact_macro_f1",
]


def _category_sort_key(cat):
    return CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else 999


def save_metrics_tsv(all_metrics: Dict, output_path: Path):
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(TSV_HEADER) + "\n")
        for model in sorted(all_metrics.keys()):
            for cat in sorted(all_metrics[model].keys(), key=_category_sort_key):
                m = all_metrics[model][cat]
                line = [cat, model, m.get("n_variants", 0), m.get("directional_accuracy", 0)]
                for d in DIRECTION_TARGETS:
                    dm = m.get(d, {})
                    line += [dm.get("sensitivity", 0), dm.get("specificity", 0), dm.get("precision", 0), dm.get("f1", 0)]

                # Exact macro F1 over the five classes
                exact_f1s = [m.get(f"exact_{c}", {}).get("f1", 0) for c in CLASSIFICATION_LABELS]
                ex_macro = sum(exact_f1s) / 5 if exact_f1s else 0
                line += [m.get("exact_accuracy", 0), f"{ex_macro:.2f}"]
                f.write("\t".join(map(str, line)) + "\n")
    print(f"Saved metrics TSV to {output_path}")
