"""
Carry over registered models whose raw results are not in this checkout,
copying their computed results from the dataset embedded in the current
dashboard/index.html.
"""

import json
from pathlib import Path

from dashboard.variant_table import NO_CALL

DASHBOARD_DIR = Path(__file__).resolve().parent


def load_previous_dataset() -> dict:
    """Dataset embedded in the current dashboard/index.html, if any."""
    html_path = DASHBOARD_DIR / "index.html"
    if not html_path.exists():
        return {}
    html = html_path.read_text(encoding="utf-8")
    marker = "window.BENCHMARK_DATA = "
    start = html.find(marker)
    if start == -1:
        return {}
    try:
        data, _ = json.JSONDecoder().raw_decode(html, start + len(marker))
    except json.JSONDecodeError:
        return {}
    return data


def carry_over_archived_models(payload: dict, missing: set[str]) -> None:
    """
    Keep registered models whose results file is not in this checkout (for
    example a run excluded by .gitignore) by copying their computed results
    from the previous build. They are flagged "archived" and lack
    variant-level comparisons against models that were added later.
    """
    previous = load_previous_dataset()
    current = {m["key"] for m in payload["models"]}
    archived = [m for m in previous.get("models", [])
                if m["key"] in missing and m["key"] not in current
                and m["key"] in previous.get("confusion_matrices", {})]

    for m in archived:
        k = m["key"]
        m = dict(m, archived=True)
        payload["models"].append(m)
        payload["confusion_matrices"][k] = previous["confusion_matrices"][k]
        for section in ("categories_data", "criteria_data", "safety_data", "cross_direction"):
            if k in previous.get(section, {}):
                payload[section][k] = previous[section][k]
        cats = payload["categories_data"].get(k, {})
        for r in cats.values():
            r.setdefault("exact_correct", round(r["exact_accuracy"] * r["count"] / 100))
            r.setdefault("directional_correct", round(r["directional_accuracy"] * r["count"] / 100))

        old_conc = previous.get("concordance", {}).get(k, {})
        payload["concordance"][k] = {o: v for o, v in old_conc.items() if o in payload["concordance"] or o == k}
        for o, v in old_conc.items():
            if o in payload["concordance"] and o != k:
                payload["concordance"][o][k] = v

        old_table = previous.get("variant_table") or {}
        old_calls = dict(zip((r[0] for r in old_table.get("rows", [])), old_table.get("calls", {}).get(k, "")))
        table = payload["variant_table"]
        table["calls"][k] = "".join(old_calls.get(r[0], NO_CALL) for r in table["rows"])
        print(f"  Carried over {m['name']} from the previous build (raw results not found).")
