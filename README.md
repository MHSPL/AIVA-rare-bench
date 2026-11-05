# AIVA Rare Bench

Benchmarks for AI variant interpretation in rare disease, with two parts:

- **Variant Classification**: ACMG/AMP 5-class calls scored against expert-panel classifications from the [ClinGen Evidence Repository (eRepo)](https://erepo.clinicalgenome.org/evrepo/). Data in `data/classification/`, runs in `results/classification/`.
- **Variant Ranking**: to be added, in `data/ranking/` and `results/ranking/`.

The dashboard (`dashboard/`, served at the repository's GitHub Pages site) shows both. The rest of this README covers Variant Classification.

## Key Results

AIVA (powered by Gemini 3 Pro) achieves **80.5% macro F1** on 8,387 expert-curated variants across 11 disease categories:

| Metric | Score |
|--------|-------|
| Overall Macro F1 | 80.5% |
| Pathogenic/Likely Pathogenic F1 | 87.4% |
| VUS F1 | 73.3% |
| Benign/Likely Benign F1 | 80.9% |

**AIVA outperforms BIAS-2015 and InterVar in 9 of 11 disease categories.**

See [BENCHMARK_FINDINGS.md](results/classification/Results_Gemini_Pro_Analysis/BENCHMARK_FINDINGS.md) for detailed analysis.

## Overview

This tool measures how well AIVA's variant classification aligns with expert panel (VCEP) curations using ACMG/AMP guidelines. It supports two primary modes of operation:

### Option 1: Standard Evaluation (Recommended)
Directly use the benchmark dataset (`data/classification/erepo_bias_intervar_merged.jsonl`): 2,894 ClinGen eRepo variants with the expert classification and the BIAS-2015 and InterVar predictions, chosen from the full 8,387-variant cohort as the variants that are hardest to classify (see [dashboard/README.md](dashboard/README.md#benchmark-variant-set)).

### Option 2: Full Reproduction
Re-generate the full 8,387-variant cohort from the BIAS-2015 supplementary files and eRepo source data (Workflow B); it is written to `data/classification/erepo_bias_intervar_full.jsonl` and does not replace the benchmark file.

### Benchmark data

| File | Contents | Use it to |
|------|----------|-----------|
| `data/classification/aiva_bench_input.jsonl` | The 2,894 variants with no classifications: `id`, `variation`, `variation_original`, `hgvs_expressions`, `variant_key`, `gene_symbol`, `disease`, `mode_of_inheritance`, `genome_build` (GRCh37, for `variant_key`), `expert_panel` | Run any model or tool on the benchmark without seeing the answers |
| `data/classification/erepo_bias_intervar_merged.jsonl` | The same 2,894 variants with the answers: eRepo expert classification and met criteria, BIAS-2015 and InterVar predictions | Score results (matched on `id`) |

To benchmark your own classifier, run it on `aiva_bench_input.jsonl` and write
one JSON line per variant with its `id`, a 5-class `classification`
(Pathogenic, Likely Pathogenic, Uncertain Significance, Likely Benign, Benign)
and, optionally, the ACMG `criteria` it applied:

```json
{"id": "d9f52cde-a2d2-4580-9f3d-af90946ce1f9", "classification": "Likely Pathogenic", "criteria": ["PM2", "PP3"]}
```

Add the file to `dashboard/models.tsv` (format `jsonl`) and rebuild; it is
scored against `erepo_bias_intervar_merged.jsonl` by `id`
(see [OPERATIONS_GUIDE.md](OPERATIONS_GUIDE.md)).

**Metrics calculated:**
- **Classification accuracy** (5-class and 3-class directional)
- **Per-class precision, recall, and F1 scores**
- **ACMG criteria agreement** (which evidence codes each tool identifies)
- **Confusion matrices** for error analysis

## Quick Start

**Requirements:** Python 3.10+ and (for the classification step only) an AIVA API token.

```bash
# 1. Setup (always)
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

#### Path A — Reproduce the analysis (no API token needed)

Regenerate all metrics and figures from the classification results bundled in
this repo. This is the fastest way to verify the pipeline works end-to-end:

```bash
python scripts/benchmark_analysis.py \
  results/classification/Results_Gemini_Pro_Analysis/classifications_google_gemini-3-pro-preview_combined.jsonl \
  --output-dir scratch/my_run
```

Outputs metrics TSVs, scorecards, and figures under `scratch/my_run/`. The
BIAS-2015 / InterVar baselines are read from the bundled
`data/classification/erepo_bias_intervar_merged.jsonl`.

#### Path B — Run AIVA classification yourself (requires your API token)

```bash
cp .env.example .env        # then edit .env and set AIVA_API_TOKEN

# Quick test on 100 variants of the 2,894-variant benchmark
python scripts/classify.py \
  --model google/gemini-3-pro-preview \
  --input data/classification/erepo_bias_intervar_merged.jsonl \
  --limit 100 \
  --output scratch/my_run/classifications.jsonl

# Then analyze your own results (no API needed for this step)
python scripts/benchmark_analysis.py scratch/my_run/classifications.jsonl --output-dir scratch/my_run
```

#### Path C — Rebuild the full cohort from source

```bash
# Download and unzip the BIAS-2015 supplementary files (Workflow B, Step 1), then:
python scripts/parse_bias_data.py --source-dir /tmp/bias_supp/bias_v2.0.0_supplementary_files
# Writes the full 8,387-variant cohort to data/classification/erepo_bias_intervar_full.jsonl;
# the 2,894-variant benchmark file is left unchanged.
```

## Setup

```bash
# Clone and enter directory
git clone https://github.com/MHSPL/AIVA-rare-bench.git
cd AIVA-rare-bench

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Configure API access
cp .env.example .env
# Edit .env with your AIVA API token
```

## Configuration

Edit `config.py` or set environment variables:

| Setting | Default | Description |
|---------|---------|-------------|
| `AIVA_API_URL` | `http://localhost:8000` | AIVA API endpoint |
| `AIVA_API_TOKEN` | - | AIVA API token for authentication |
| `DEFAULT_MODEL` | `google/gemini-3-pro-preview` | LLM model for classification |
| `DEFAULT_ASSEMBLY` | `GRCh37` | Genome assembly (eRepo baseline is GRCh37/hg19) |
| `CONCURRENCY_LIMIT` | `5` | Max concurrent API requests (respects AIVA rate limit) |

## Project Structure

```
AIVA-rare-bench/
├── config.py                    # Configuration and paths
├── requirements.txt             # Python dependencies
├── .env.example                 # Environment template
├── OPERATIONS_GUIDE.md          # Adding or updating a model (run, register, rebuild, PR)
├── index.html                   # GitHub Pages entry; forwards to dashboard/
├── dashboard/                   # AIVA Rare Bench dashboard (see dashboard/README.md)
│   ├── models.tsv               # Registry: which results files appear, and how
│   ├── build_data.py            # Scores every registered model -> data/benchmark_data.json
│   ├── generate_html.py         # template.html + data -> index.html (static page)
│   └── index.html               # The built dashboard
├── data/
│   ├── classification/
│   │   ├── aiva_bench_input.jsonl            # Input: 2,894 variants, no classifications
│   │   └── erepo_bias_intervar_merged.jsonl  # Answers: expert class + BIAS-2015 + InterVar
│   └── ranking/                 # Variant ranking inputs and answers (to come)
├── results/
│   ├── classification/          # One folder per model run
│   │   ├── Results_<Model>_Analysis/  # AIVA runs: <model>_full.jsonl, metrics, figures
│   │   └── Claude_code_results/       # Claude Code (agent) classifications
│   └── ranking/                 # Variant ranking runs (to come)
├── scratch/                     # Ad-hoc runs and analyses (git-ignored)
└── scripts/
    ├── classify.py              # Run a model through the AIVA API (resumes)
    ├── retry_errors.py          # Re-run failed classifications
    ├── ledger.py                # Coverage: what is classified, what is left
    ├── update_dashboard.py      # Rebuild the dashboard data and page
    ├── sample_stratified.py     # Stratified sample for a pilot run
    ├── parse_bias_data.py       # Rebuild the full cohort from raw source
    ├── benchmark_analysis.py    # End-to-end analysis pipeline (metrics + figures)
    ├── generate_results.py      # Metrics calculation & figure generation
    ├── metrics.py               # Shared metric utilities
    └── review_results.R         # Step through a run by hand in R
```

## Usage

### Workflow A: Run a Model and Add It to the Dashboard

Requires an AIVA API token in `.env`. `classify.py` resumes, so re-running
the same command picks up where it stopped; use `--limit` for a quick test.

```bash
python scripts/classify.py \
  --model google/gemini-3-pro-preview \
  --input data/classification/erepo_bias_intervar_merged.jsonl \
  --output results/classification/Results_<Model>_Analysis/<model>_full.jsonl \
  --concurrency 5
```

Then add the results file to `dashboard/models.tsv` and run
`python scripts/update_dashboard.py`; the dashboard computes F1, confusion
matrices, safety and agreement for every registered model. The full workflow,
including the retry loop and the pull request, is in
[OPERATIONS_GUIDE.md](OPERATIONS_GUIDE.md).

---

### Workflow B: 4-Way Tool Comparison (AIVA vs InterVar vs BIAS-2015)

The benchmark file already holds the BIAS-2015 and InterVar predictions for its 2,894 variants. These steps rebuild the full 8,387-variant cohort from the original sources, only needed to audit or re-derive it.

#### Step 1: Download BIAS-2015 Supplementary Data

Download pre-computed predictions from the BIAS-2015 paper:

**Option A: Direct download (no AWS CLI required)**

Download: https://bias-2015.s3.amazonaws.com/bias_v2.0.0_supplementary_files.zip

**Option B: AWS CLI**

```bash
aws s3 cp s3://bias-2015/bias_v2.0.0_supplementary_files.zip /tmp/
```

Then extract and copy files:

```bash
unzip /tmp/bias_v2.0.0_supplementary_files.zip -d /tmp/bias_supp/

# Copy required files to project
cp /tmp/bias_supp/bias_v2.0.0_supplementary_files/erepo_feb_2nd.tabbed data/classification/bias_2015/
cp /tmp/bias_supp/bias_v2.0.0_supplementary_files/erepo_feb_2nd.hg19_multianno.txt.intervar data/classification/bias_2015/
cp /tmp/bias_supp/bias_v2.0.0_supplementary_files/erepo_feb_2nd_bitscopic_classifications.tsv data/classification/bias_2015/
cp /tmp/bias_supp/bias_v2.0.0_supplementary_files/erepo_feb_2nd_sorted.vcf data/classification/bias_2015/
```

**Files needed:**
| File | Description |
|------|-------------|
| `erepo_feb_2nd.tabbed` | eRepo ground truth (expert classifications) |
| `erepo_feb_2nd.hg19_multianno.txt.intervar` | InterVar predictions |
| `erepo_feb_2nd_bitscopic_classifications.tsv` | BIAS-2015 predictions |
| `erepo_feb_2nd_sorted.vcf` | Authoritative variant list (VCF format) |

#### Step 2: Parse and Unify Data

```bash
python scripts/parse_bias_data.py
```

This creates `data/classification/erepo_bias_intervar_full.jsonl` containing 8,387 variants with:
- Ground truth classifications from eRepo
- InterVar predictions
- BIAS-2015 predictions
- All in AIVA-compatible format (hg19/GRCh37)

#### Step 3: Run AIVA and Compare

Run AIVA on the benchmark file as in Workflow A, then rebuild the dashboard.
BIAS-2015 and InterVar are registered in `dashboard/models.tsv`, so the
leaderboard scores AIVA, the Claude Code agents and both rule-based tools on
the same variants.

## Script Reference

| Script | Purpose |
|--------|---------|
| `classify.py` | Run a model through the AIVA API; resumes an interrupted run |
| `retry_errors.py` | Re-run failed classifications in a results file |
| `ledger.py` | Report what is classified so far, and what is left |
| `update_dashboard.py` | Rebuild the dashboard data and page |
| `sample_stratified.py` | Draw a stratified sample for a pilot run |
| `parse_bias_data.py` | Parse raw eRepo/BIAS-2015/InterVar files, create the unified dataset |
| `benchmark_analysis.py` | End-to-end analysis pipeline with figures |
| `generate_results.py` | Metrics calculation and figure generation |
| `metrics.py` | Shared metric utilities |
| `review_results.R` | Step through a run by hand in R |

## Metrics Explained

### Classification Accuracy

| Metric | Description |
|--------|-------------|
| **Exact Match (5-class)** | % where AIVA matches expert classification exactly |
| **Directional (3-class)** | % correct when collapsed to Pathogenic/VUS/Benign |

### Per-Class Metrics

For each of the 5 ACMG classifications:
- **Precision**: When AIVA predicts this class, how often is it correct?
- **Recall**: Of all variants in this class, how many does AIVA identify?
- **F1**: Harmonic mean of precision and recall

### Criteria Analysis

| Metric | Description |
|--------|-------------|
| **Avg Jaccard (normalized)** | Overlap of ACMG criteria (after normalizing modifiers like `_Supporting`) |
| **Direction Agreement** | % where AIVA's evidence leans the same direction (pathogenic vs benign) |
| **Per-Criterion P/R/F1** | How well AIVA identifies each specific criterion (PM2, PP3, etc.) |

### Evidence Direction

Criteria are grouped by evidence direction:
- **Pathogenic**: PVS1, PS1-4, PM1-6, PP1-5
- **Benign**: BA1, BS1-4, BP1-7

## Rate Limiting

The classify script handles rate limits automatically:
- Exponential backoff on errors
- Special handling for 429 (rate limit) with longer delays
- Configurable concurrency to avoid overwhelming the API

If you encounter rate limits, reduce concurrency:
```bash
python scripts/classify.py --concurrency 1 ...
```

## License

MIT License - see LICENSE file for details.

## Contact

Arun Boddapati, MHSPL Research Team
