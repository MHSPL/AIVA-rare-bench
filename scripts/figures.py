"""
Classification figures for RareDx Bench: confusion matrices and the
misclassification breakdown (figure 10).
"""

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.categories import get_short_name
from scripts.metrics import get_direction


def plot_confusion_matrix(cm, classes, title, output_path):
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    plt.figure(figsize=(10, 8))
    plt.imshow(cm_norm, interpolation='nearest', cmap=plt.cm.Blues)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, ha="right")
    plt.yticks(tick_marks, classes)

    thresh = cm_norm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, f"{cm[i, j]}\n({cm_norm[i, j]:.1%})", ha="center", color="white" if cm_norm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('Ground Truth')
    plt.xlabel('Predicted')
    plt.savefig(output_path, dpi=300)
    plt.close()


def _best_model(aiva_data: dict):
    """The model with the highest directional accuracy."""
    best_model = None
    best_acc = 0
    for model, variants in aiva_data.items():
        correct = sum(1 for v in variants
                      if get_direction(v.get("ground_truth", "")) == get_direction(v.get("aiva_classification", ""))
                      and not v.get("error"))
        total = sum(1 for v in variants if not v.get("error"))
        acc = correct / total if total > 0 else 0
        if acc > best_acc:
            best_acc = acc
            best_model = model
    return best_model


# Internal keys match get_direction() output; labels are for display.
DIR_KEYS = ['pathogenic', 'vus', 'benign']
DISPLAY_LABELS = {'pathogenic': 'P/LP', 'vus': 'VUS', 'benign': 'B/LB'}
# Colors matched to the Sankey flow diagram palette
DIR_COLORS = {'pathogenic': '#c44e52', 'vus': '#edc948', 'benign': '#4c72b0'}
# Each misclassification type gets a shade of its source (true) class color
DONUT_COLOR_MAP = {
    ('benign', 'vus'):        '#6a9fd8',  # lighter blue (B/LB source)
    ('benign', 'pathogenic'): '#8bb8e8',  # lightest blue (B/LB source, rare)
    ('vus', 'pathogenic'):    '#f5de78',  # lighter gold (VUS source)
    ('vus', 'benign'):        '#f0d04a',  # darker gold (VUS source)
    ('pathogenic', 'vus'):    '#d47a7d',  # lighter red (P/LP source)
    ('pathogenic', 'benign'): '#e0a3a5',  # lightest red (P/LP source, rare)
}


def _plot_distribution(ax, correct_types, misclass_types, model_name):
    """Stacked bars: for each true class, where its predictions went."""
    x = np.arange(len(DIR_KEYS))
    bottom = np.zeros(len(DIR_KEYS))

    for pred_class in DIR_KEYS:
        heights = []
        for true_class in DIR_KEYS:
            if true_class == pred_class:
                heights.append(correct_types.get(true_class, 0))
            else:
                heights.append(misclass_types.get((true_class, pred_class), 0))

        ax.bar(x, heights, bottom=bottom, label=f'Pred: {DISPLAY_LABELS[pred_class]}',
               color=DIR_COLORS[pred_class], edgecolor='black', linewidth=0.5, alpha=0.8)
        bottom += heights

    ax.set_xticks(x)
    ax.set_xticklabels([f'True:\n{DISPLAY_LABELS[l]}' for l in DIR_KEYS])
    ax.set_ylabel('Number of Variants')
    ax.set_title(f'Classification Distribution by True Class\n({model_name})', fontweight='bold')
    ax.legend(loc='upper right')


def _plot_error_breakdown(ax, misclass_types, model_name):
    """Pie chart of misclassification types, keyed by a legend."""
    error_counts, error_labels, error_colors = [], [], []
    for (true_c, pred_c), count in sorted(misclass_types.items(), key=lambda x: -x[1]):
        if count > 0:
            error_counts.append(count)
            error_labels.append(f'{DISPLAY_LABELS[true_c]} → {DISPLAY_LABELS[pred_c]} (n={count})')
            error_colors.append(DONUT_COLOR_MAP.get((true_c, pred_c), '#cccccc'))

    if error_counts:
        # A legend avoids overlapping labels; wedges stay flush so the center is filled.
        wedges, texts, autotexts = ax.pie(error_counts, autopct='%1.1f%%',
                                          colors=error_colors, explode=[0.0] * len(error_counts),
                                          startangle=140, pctdistance=0.85)
        ax.legend(wedges, error_labels, title="Error Types", loc="center left",
                  bbox_to_anchor=(1, 0, 0.5, 1), fontsize=8)
        ax.set_title(f'Misclassification Breakdown\n({model_name})', fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No errors!', ha='center', va='center', fontsize=14)
        ax.set_title('Misclassification Breakdown', fontweight='bold')


def plot_error_analysis(aiva_data: dict, output_dir: Path):
    """Analyze and visualize misclassification patterns."""
    best_model = _best_model(aiva_data)
    if not best_model:
        print("No AIVA model found for error analysis")
        return

    # Collect misclassifications for the best AIVA model
    misclass_types = defaultdict(int)
    correct_types = defaultdict(int)
    for v in aiva_data[best_model]:
        gt = get_direction(v.get("ground_truth", ""))
        pred = get_direction(v.get("aiva_classification", ""))
        if gt == "unknown" or pred == "unknown" or v.get("error"):
            continue
        if gt == pred:
            correct_types[gt] += 1
        else:
            misclass_types[(gt, pred)] += 1

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    model_name = get_short_name(best_model)
    _plot_distribution(axes[0], correct_types, misclass_types, model_name)
    _plot_error_breakdown(axes[1], misclass_types, model_name)
    plt.tight_layout()

    output_path = output_dir / 'figure10_error_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure10_error_analysis.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
