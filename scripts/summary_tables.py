"""
Scorecard tables for RareDx Bench, drawn with matplotlib from
metrics_by_category.tsv: the per-category summary and the combined
(overall plus category) summary.
"""

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

TEXT_COLOR = "#202124"
MUTED_COLOR = "#5f6368"
AIVA_BG_COLOR = "#dbeafe"
OTHERS_BG_COLOR = "#f1f3f4"
LINE_COLOR = "#e0e0e0"
SEPARATOR_COLOR = "#9aa0a6"


def macro_f1(row) -> float:
    """Macro F1 as the mean of the P/LP, VUS and B/LB F1 scores."""
    return (row['plp_f1'] + row['vus_f1'] + row['blb_f1']) / 3


def _model_rows(df, key):
    """Rows for one model: substring match for AIVA (gemini) runs, exact otherwise."""
    return df[df['model'].str.contains(key)] if 'gemini' in key else df[df['model'] == key]


def category_macro_f1_table(df) -> pd.DataFrame:
    """Macro F1 (%) per category for Gemini 3 Pro, BIAS-2015 and InterVar."""
    cats = df['category'].unique()
    data = {"Category": cats}
    model_map = {'gemini-3-pro-preview': 'Gemini 3 Pro', 'BIAS-2015': 'BIAS-2015', 'InterVar': 'InterVar'}
    for m_key, m_disp in model_map.items():
        scores = []
        for c in cats:
            row = _model_rows(df[df['category'] == c], m_key)
            scores.append(macro_f1(row.iloc[0]) * 100 if not row.empty else 0)
        data[m_disp] = scores
    return pd.DataFrame(data)


def _column_positions(model_columns):
    """Normalized column widths and left edges: a wide category column, then models."""
    cat_width = 0.40
    model_width = 0.60 / len(model_columns) if model_columns else 0.1
    col_widths = [cat_width] + [model_width] * len(model_columns)
    total_width = sum(col_widths)
    norm_widths = [w / total_width for w in col_widths]
    x_pos = [0]
    for w in norm_widths:
        x_pos.append(x_pos[-1] + w)
    return norm_widths, x_pos


def generate_category_summary_table(summary_df, output_dir):
    """Category x model table: AIVA columns on blue, others on grey."""
    summary_df = summary_df.sort_values('Category')
    model_columns = [col for col in summary_df.columns if col != 'Category']
    aiva_models = [col for col in model_columns if 'Gemini' in col]
    other_models = [col for col in model_columns if col not in aiva_models]

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    header_height = 0.08
    row_height = 0.055
    norm_widths, x_pos = _column_positions(model_columns)
    table_top = 0.82
    table_bottom = table_top - (len(summary_df) + 1) * row_height
    header_top = table_top + 0.005

    # AIVA background and header
    aiva_indices = [i + 1 for i, m in enumerate(model_columns) if m in aiva_models]
    if aiva_indices:
        start, end = min(aiva_indices), max(aiva_indices)
        bg_x, bg_w = x_pos[start], sum(norm_widths[start: end + 1])
        ax.add_patch(plt.Rectangle((bg_x, table_bottom), bg_w, (table_top - table_bottom) + row_height * 0.5, color=AIVA_BG_COLOR, transform=ax.transAxes, zorder=0))
        ax.add_patch(plt.Rectangle((bg_x, header_top), bg_w, header_height * 0.7, color=AIVA_BG_COLOR, transform=ax.transAxes, zorder=0))
        ax.text(bg_x + bg_w / 2, header_top + header_height * 0.35, "AIVA", ha='center', va='center', fontweight='bold', fontsize=13, color="#1967d2", transform=ax.transAxes)

    # Others background, header and separator
    other_indices = [i + 1 for i, m in enumerate(model_columns) if m in other_models]
    if other_indices:
        start, end = min(other_indices), max(other_indices)
        bg_x, bg_w = x_pos[start], sum(norm_widths[start: end + 1])
        ax.add_patch(plt.Rectangle((bg_x, header_top), bg_w, header_height * 0.7, color=OTHERS_BG_COLOR, transform=ax.transAxes, zorder=0))
        ax.text(bg_x + bg_w / 2, header_top + header_height * 0.35, "Others", ha='center', va='center', fontweight='bold', fontsize=13, color=MUTED_COLOR, transform=ax.transAxes)
        if aiva_indices:
            sep_x = x_pos[start]
            ax.plot([sep_x, sep_x], [table_bottom, header_top + header_height * 0.7], color=SEPARATOR_COLOR, linewidth=1.5, transform=ax.transAxes, zorder=2)

    # Column headers
    headers = ['Category'] + model_columns
    current_y = table_top - row_height
    ax.plot([0, 1], [current_y, current_y], color='#333', linewidth=1.0, transform=ax.transAxes)
    for i, h in enumerate(headers):
        align = 'left' if i == 0 else 'center'
        x_txt = x_pos[i] + 0.02 if align == 'left' else x_pos[i] + norm_widths[i] / 2
        ax.text(x_txt, current_y + row_height / 2, h, transform=ax.transAxes, fontweight='bold', fontsize=11, ha=align, va='center', color=TEXT_COLOR)

    # Rows: the best score in each row is bold
    for idx, row in summary_df.iterrows():
        current_y -= row_height
        ax.plot([0, 1], [current_y, current_y], color=LINE_COLOR, linewidth=0.5, transform=ax.transAxes)
        ax.text(x_pos[0] + 0.02, current_y + row_height / 2, row['Category'], transform=ax.transAxes, fontsize=10, fontweight='medium', va='center', color=TEXT_COLOR)
        max_val = max(row[m] for m in model_columns)
        for i, m in enumerate(model_columns):
            val = row[m]
            weight, size = ('bold', 11) if val == max_val else ('normal', 10)
            x_txt = x_pos[i + 1] + norm_widths[i + 1] / 2
            ax.text(x_txt, current_y + row_height / 2, f"{val:.1f}%", transform=ax.transAxes, ha='center', va='center', fontweight=weight, fontsize=size, color=TEXT_COLOR)

    ax.text(0, table_bottom - row_height * 0.8, "Source: Internal Benchmark.", fontsize=8, color=MUTED_COLOR, transform=ax.transAxes)

    out = Path(output_dir)
    plt.savefig(out / "figure_category_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(out / "figure_category_summary.pdf", bbox_inches='tight')
    print(f"Saved category summary to {out}/figure_category_summary.png")


COMBINED_MODELS = {'gemini-3-pro-preview': 'AIVA*', 'BIAS-2015': 'BIAS-2015', 'InterVar': 'InterVar'}


def build_combined_table(df) -> pd.DataFrame:
    """Overall rows (variant weighted F1s), a spacer, then one Macro F1 row per category."""
    cat_data = {"Category": []}
    cat_n_variants = {}
    for disp in COMBINED_MODELS.values():
        cat_data[disp] = []

    for cat in sorted(df['category'].unique()):
        cat_df = df[df['category'] == cat]
        gemini_row = cat_df[cat_df['model'].str.contains('gemini')]
        n_var = int(gemini_row['n_variants'].values[0]) if not gemini_row.empty else 0
        cat_n_variants[cat] = n_var
        cat_data["Category"].append(f"{cat} (n={n_var})")
        for key, disp in COMBINED_MODELS.items():
            row = _model_rows(cat_df, key)
            cat_data[disp].append(macro_f1(row.iloc[0]) * 100 if not row.empty else 0)

    total_variants = sum(cat_n_variants.values())
    overall_data = {"Category": [
        f"Overall (n={total_variants})",
        "Pathogenic / Likely Pathogenic",
        "Variant of Uncertain Significance",
        "Benign / Likely Benign",
        ""  # Spacer
    ]}
    for key, disp in COMBINED_MODELS.items():
        model_df = _model_rows(df, key)
        if model_df.empty:
            overall_data[disp] = [0, 0, 0, 0, '']
            continue
        total_n = model_df['n_variants'].sum()
        w_plp = (model_df['plp_f1'] * model_df['n_variants']).sum() / total_n
        w_vus = (model_df['vus_f1'] * model_df['n_variants']).sum() / total_n
        w_blb = (model_df['blb_f1'] * model_df['n_variants']).sum() / total_n
        w_macro = (w_plp + w_vus + w_blb) / 3
        overall_data[disp] = [w_macro * 100, w_plp * 100, w_vus * 100, w_blb * 100, '']

    return pd.concat([pd.DataFrame(overall_data), pd.DataFrame(cat_data)], ignore_index=True)


def generate_combined_summary(metrics_tsv_path, output_dir):
    """Generate combined summary figure (Overall + Category) using exact TSV values."""
    summary_df = build_combined_table(pd.read_csv(metrics_tsv_path, sep='\t'))
    model_columns = [col for col in summary_df.columns if col != 'Category']
    aiva_models = ['AIVA*']
    other_models = [col for col in model_columns if col not in aiva_models]

    fig, ax = plt.subplots(figsize=(12, 14))
    ax.axis('off')

    header_height = 0.06
    row_height = 0.05
    norm_w, x_pos = _column_positions(model_columns)
    table_top = 0.88
    # Extend the bottom slightly so the shading covers the last row fully
    table_bottom = table_top - (len(summary_df) + 1.5) * row_height
    top_of_header = table_top + 0.005 + header_height * 0.7

    aiva_indices = [i + 1 for i, m in enumerate(model_columns) if m in aiva_models]
    if aiva_indices:
        s, e = min(aiva_indices), max(aiva_indices)
        bg_x, bg_w = x_pos[s], sum(norm_w[s:e + 1])
        rect_h = (table_top - table_bottom) + header_height
        ax.add_patch(plt.Rectangle((bg_x, table_bottom), bg_w, rect_h, color=AIVA_BG_COLOR, transform=ax.transAxes, zorder=0))
        # Drawn twice in the original layout; kept so the figure stays pixel identical.
        for _ in range(2):
            ax.add_patch(plt.Rectangle((bg_x, table_bottom), bg_w, top_of_header - table_bottom, color=AIVA_BG_COLOR, transform=ax.transAxes, zorder=0))

    # Separator between AIVA and the other tools (no background for the others)
    other_indices = [i + 1 for i, m in enumerate(model_columns) if m in other_models]
    if other_indices and aiva_indices:
        sep_x = x_pos[min(other_indices)]
        ax.plot([sep_x, sep_x], [table_bottom, top_of_header], color=SEPARATOR_COLOR, linewidth=1.5, transform=ax.transAxes, zorder=2)

    # Column headers, AIVA in blue
    headers = ['Benchmark (F1 Score) / Category'] + model_columns
    current_y = table_top - row_height
    ax.plot([0, 1], [current_y, current_y], color='#333', linewidth=1, transform=ax.transAxes)
    for i, h in enumerate(headers):
        align = 'left' if i == 0 else 'center'
        x = x_pos[i] + (0.02 if align == 'left' else norm_w[i] / 2)
        h_color = "#1967d2" if "AIVA" in h else TEXT_COLOR
        ax.text(x, current_y + row_height / 2, h, transform=ax.transAxes, fontweight='bold', fontsize=11, ha=align, va='center', color=h_color)

    for idx, row in summary_df.iterrows():
        current_y -= row_height
        if row['Category'] == "":
            continue  # Spacer
        ax.plot([0, 1], [current_y, current_y], color=LINE_COLOR, linewidth=0.5, transform=ax.transAxes)
        _draw_row_label(ax, row['Category'], idx, x_pos, current_y + row_height / 2)

        vals = [row[m] for m in model_columns]
        if all(isinstance(v, (int, float)) for v in vals):
            max_v = max(vals)
            for i, m in enumerate(model_columns):
                val = row[m]
                w, s = ('bold', 11) if val == max_v else ('normal', 10)
                ax.text(x_pos[i + 1] + norm_w[i + 1] / 2, current_y + row_height / 2, f"{val:.1f}%", transform=ax.transAxes, ha='center', va='center', fontweight=w, fontsize=s, color=TEXT_COLOR)

    # "Category Breakdown" sits in the second spacer row (rows 0 to 3 are overall, 4 is a spacer)
    cat_header_y = table_top - (6 * row_height)
    ax.text(x_pos[0], cat_header_y + row_height * 0.5, "Category Breakdown", transform=ax.transAxes, fontsize=11, fontweight='bold', color='#444', va='center')

    eval_date = datetime.now().strftime("%b %d, %Y")
    ax.text(0, table_bottom - row_height * 0.8, f"Source: Internal Benchmark. Evaluation date: {eval_date}.", fontsize=8, color=MUTED_COLOR, transform=ax.transAxes)
    ax.text(0, table_bottom - row_height * 1.5, "* AIVA powered by Gemini 3 Pro", fontsize=7, color=MUTED_COLOR, transform=ax.transAxes, fontstyle='italic')

    out = Path(output_dir)
    plt.savefig(out / "figure_combined_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(out / "figure_combined_summary.pdf", bbox_inches='tight')
    print(f"Saved combined summary to {out}/figure_combined_summary.png")


def _draw_row_label(ax, cat_text, idx, x_pos, y):
    """Row label: bold for benchmark and category rows, with any (n=...) count right aligned."""
    is_bench = idx < 4
    is_category = idx >= 5  # After the spacer row
    weight = 'bold' if (is_bench or is_category) else 'medium'
    label = dict(transform=ax.transAxes, fontsize=10, fontweight=weight, va='center', color=TEXT_COLOR)

    if "(n=" not in cat_text:
        ax.text(x_pos[0] + 0.02, y, cat_text, **label)
        return
    main_text, count = cat_text.split(" (n=")[:2]
    # The label is drawn twice in the original layout; kept so the figure stays pixel identical.
    ax.text(x_pos[0] + 0.02, y, main_text, **label)
    cat_col_w = x_pos[1] - x_pos[0]
    ax.text(x_pos[0] + cat_col_w - 0.05, y, "(n=" + count, transform=ax.transAxes, fontsize=9, fontstyle='italic', fontweight='normal', ha='right', va='center', color=MUTED_COLOR)
    ax.text(x_pos[0] + 0.02, y, main_text, **label)
