"""
Gene prioritization benchmark (simulated patients) for the benchmark analysis.

Scores predicted genes against each patient's true and distractor genes,
splitting hits into seed (knowledge graph or candidate) and novel predictions,
and writes gene_prioritization_metrics.json and gene_prioritization_summary.tsv.
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.metrics import calculate_precision_recall_f1


def normalize_gene_symbol(gene: str) -> str:
    if not gene:
        return ""
    return gene.strip().upper()


def _parse_benchmark_entry(entry: dict) -> tuple[str, dict]:
    """One simulated patient: (patient_id, gene sets and categories)."""
    gene_id_to_symbol = {}
    for v in entry.get("variants", []):
        gene_id = v.get("gene_id")
        if gene_id:
            gene_id_to_symbol[str(gene_id)] = v.get("gene_symbol") or gene_id

    def to_symbol_set(gene_ids):
        return {normalize_gene_symbol(gene_id_to_symbol.get(str(gid), str(gid))) for gid in gene_ids or []}

    true_genes = to_symbol_set(entry.get("true_genes", []))
    distractor_genes = to_symbol_set(entry.get("distractor_genes", {}).keys())
    all_candidates = entry.get("all_candidate_genes", [])
    candidate_genes = to_symbol_set(all_candidates) if all_candidates else (true_genes | distractor_genes)

    return str(entry.get("id")), {
        "true_genes": true_genes,
        "distractor_genes": distractor_genes,
        "candidate_genes": candidate_genes,
        "category": entry.get("category", "unknown"),
        "broad_category": entry.get("broad_category", "unknown"),
    }


def load_gene_benchmark(benchmark_path: Path) -> dict:
    """Load simulated benchmark JSONL into a patient_id -> benchmark dict."""
    benchmark = {}
    with open(benchmark_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                patient_id, data = _parse_benchmark_entry(json.loads(line))
                benchmark[patient_id] = data
    return benchmark


def load_gene_benchmark_from_index(db_path: Path) -> dict:
    """Load benchmark entries from a SQLite index created by index_benchmark.py"""
    import sqlite3
    if not db_path.exists():
        raise FileNotFoundError(f"Index DB not found: {db_path}")
    bench = {}
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("SELECT id, data FROM patients")
    for row in cur:
        try:
            obj = json.loads(row[1])
        except Exception:
            continue
        patient_id, data = _parse_benchmark_entry(obj)
        bench[patient_id] = data
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
            predicted = entry.get("predicted_genes", []) or []
            kg_evidence = entry.get("kg_evidence") or entry.get("kg_evidence_json") or {}
            # kg_evidence_json may be a JSON string
            if isinstance(kg_evidence, str) and kg_evidence.strip():
                try:
                    kg_evidence = json.loads(kg_evidence)
                except Exception:
                    kg_evidence = {}
            predictions[str(entry.get("patient_id"))] = {
                "predicted": [normalize_gene_symbol(g) for g in predicted if g],
                "kg_evidence": kg_evidence
            }
    return predictions


SUMMED_FIELDS = ["tp", "fp", "fn", "fp_in_candidates", "fp_outside_candidates"]
OPTIONAL_FIELDS = ["seed_hits", "novel_hits", "seed_predictions", "novel_predictions"]


def summarize_gene_metrics(rows: list[dict]) -> dict:
    totals = {"n_patients": 0, **{k: 0 for k in SUMMED_FIELDS}, "hits": 0, "predicted_genes_total": 0,
              **{k: 0 for k in OPTIONAL_FIELDS}}

    for r in rows:
        totals["n_patients"] += 1
        for k in SUMMED_FIELDS:
            totals[k] += r[k]
        totals["hits"] += 1 if r["tp"] > 0 else 0
        totals["predicted_genes_total"] += r["predicted_count"]
        for k in OPTIONAL_FIELDS:
            totals[k] += r.get(k, 0)

    n = totals["n_patients"]
    precision, recall, f1 = calculate_precision_recall_f1(totals["tp"], totals["fp"], totals["fn"])
    totals["precision"] = round(precision, 4)
    totals["recall"] = round(recall, 4)
    totals["f1"] = round(f1, 4)
    totals["hit_rate"] = round(totals["hits"] / n, 4) if n else 0
    totals["avg_predicted_genes"] = round(totals["predicted_genes_total"] / n, 2) if n else 0
    totals["seed_hit_rate"] = round(totals["seed_hits"] / totals["seed_predictions"], 4) if totals["seed_predictions"] else None
    totals["novel_hit_rate"] = round(totals["novel_hits"] / totals["novel_predictions"], 4) if totals["novel_predictions"] else None
    return totals


def _seed_genes(kg_evidence, candidate_set) -> set:
    """Genes the knowledge graph marks as seeds; the candidates when it marks none."""
    seed_set = set()
    if isinstance(kg_evidence, dict):
        for g, ev in kg_evidence.items():
            status = ""
            if isinstance(ev, dict):
                status = str(ev.get("status", "")).lower()
            elif isinstance(ev, str):
                status = ev.lower()
            if "seed" in status:
                seed_set.add(normalize_gene_symbol(g))
    return seed_set or {normalize_gene_symbol(g) for g in candidate_set}


def score_patient(patient_id: str, b: dict, pred_entry) -> dict:
    """TP/FP/FN and seed versus novel hits for one patient."""
    predicted = pred_entry.get("predicted", []) if isinstance(pred_entry, dict) else pred_entry or []
    kg_evidence = pred_entry.get("kg_evidence", {}) if isinstance(pred_entry, dict) else {}

    pred_set = set(predicted)
    true_set = b["true_genes"]
    candidate_set = b["candidate_genes"]
    seed_set = _seed_genes(kg_evidence, candidate_set)

    return {
        "patient_id": patient_id,
        "category": b["category"],
        "broad_category": b["broad_category"],
        "tp": len(pred_set & true_set),
        "fp": len(pred_set - true_set),
        "fn": len(true_set - pred_set),
        "fp_in_candidates": len(pred_set & (candidate_set - true_set)),
        "fp_outside_candidates": len(pred_set - candidate_set),
        "predicted_count": len(pred_set),
        "seed_predictions": len(pred_set & seed_set),
        "novel_predictions": len(pred_set - seed_set),
        "seed_hits": len((pred_set & seed_set) & true_set),
        "novel_hits": len((pred_set - seed_set) & true_set),
    }


SUMMARY_COLUMNS = ["n_patients", "precision", "recall", "f1", "hit_rate", "avg_predicted_genes",
                   "fp_in_candidates", "fp_outside_candidates"]


def write_gene_outputs(overall, category_summary, broad_summary, output_dir: Path):
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
        writer.writerow(["group_type", "group_name"] + SUMMARY_COLUMNS)
        writer.writerow(["overall", "overall"] + [overall[c] for c in SUMMARY_COLUMNS])
        for group_type, summaries in (("category", category_summary), ("broad_category", broad_summary)):
            for name, summary in sorted(summaries.items()):
                writer.writerow([group_type, name] + [summary[c] for c in SUMMARY_COLUMNS])

    print(f"  [8/8] Wrote gene prioritization metrics to {metrics_path.name}")
    print(f"        Wrote gene prioritization summary to {summary_path.name}")


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
        row = score_patient(patient_id, b, predictions.get(patient_id, {}))
        if row["predicted_count"] == 0:
            missing_predictions += 1
        rows.append(row)

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
    write_gene_outputs(overall, category_summary, broad_summary, output_dir)
