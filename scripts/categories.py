"""
Disease categories for RareDx Bench.

Maps each ClinGen Expert Panel (VCEP) to a broader disease category, and maps
model identifiers to short display names. Single source of truth for the
analysis scripts and the dashboard build.
"""

TOOL_SHORT_NAMES = {
    'anthropic/claude-opus-4-5-20251101': 'AIVA-Claude',
    'google/gemini-3-pro-preview': 'AIVA-Gemini',
    'openai/gpt-5.2': 'AIVA-GPT',
    'claude-opus-4-5-20251101': 'AIVA-Claude',
    'gemini-3-pro-preview': 'AIVA-Gemini',
    'gpt-5.2': 'AIVA-GPT',
    'BIAS-2015': 'BIAS-2015',
    'InterVar': 'InterVar',
}


def get_short_name(tool_name: str) -> str:
    """Get short display name for a tool."""
    name_str = str(tool_name)
    for key, short in TOOL_SHORT_NAMES.items():
        if key in name_str or name_str in key:
            return short
    return name_str.split('/')[-1][:15]


VCEP_TO_CATEGORY = {
    # Cancer
    "Myeloid Malignancy VCEP": "Cancer Predisposition", "PTPN11 VCEP": "RASopathies",
    "CDH1 VCEP": "Cancer Predisposition", "PTEN VCEP": "Cancer Predisposition",
    "TP53 VCEP": "Cancer Predisposition", "VHL VCEP": "Cancer Predisposition",
    "Hereditary Breast, Ovarian and Pancreatic Cancer VCEP": "Cancer Predisposition",
    "ENIGMA BRCA1 and BRCA2 VCEP": "Cancer Predisposition",
    "InSiGHT Hereditary Colorectal Cancer/Polyposis VCEP": "Cancer Predisposition",
    "DICER1 and miRNA-Processing Gene VCEP": "Cancer Predisposition",
    # Metabolic
    "Monogenic Diabetes VCEP": "Metabolic Disorders", "Phenylketonuria VCEP": "Metabolic Disorders",
    "Familial Hypercholesterolemia VCEP": "Metabolic Disorders", "Lysosomal Diseases VCEP": "Metabolic Disorders",
    "Cerebral Creatine Deficiency Syndromes VCEP": "Metabolic Disorders", "ACADVL VCEP": "Metabolic Disorders",
    "Mitochondrial Diseases VCEP": "Metabolic Disorders",
    # Cardio
    "Cardiomyopathy VCEP": "Cardiovascular Disorders", "FBN1 VCEP": "Cardiovascular Disorders",
    "Pulmonary Hypertension VCEP": "Cardiovascular Disorders",
    "Potassium Channel Arrhythmia VCEP": "Cardiovascular Disorders",
    # Heme
    "Platelet Disorders VCEP": "Hematologic Disorders", "von Willebrand Disease  VCEP": "Hematologic Disorders",
    "Coagulation Factor Deficiency VCEP": "Hematologic Disorders", "Thrombosis VCEP": "Hematologic Disorders",
    # Immuno
    "Severe Combined Immunodeficiency Disease  VCEP": "Immunodeficiency Disorders",
    "Antibody Deficiencies VCEP": "Immunodeficiency Disorders",
    # Neuro
    "Rett and Angelman-like Disorders VCEP": "Neurodevelopmental Disorders",
    # Neuromuscular
    "Malignant Hyperthermia Susceptibility VCEP": "Neuromuscular Disorders",
    "Congenital Myopathies VCEP": "Neuromuscular Disorders",
    "Limb Girdle Muscular Dystrophy VCEP": "Neuromuscular Disorders",
    # Optho
    "Glaucoma VCEP": "Ophthalmologic Disorders",
    "Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP": "Ophthalmologic Disorders",
    "X-linked Inherited Retinal Disease VCEP": "Ophthalmologic Disorders", "ABCA4 VCEP": "Ophthalmologic Disorders",
    # Hearing
    "Hearing Loss VCEP": "Hearing Disorders",
    # RAS
    "RASopathy VCEP": "RASopathies", "Brain Malformations VCEP": "RASopathies",
    # Vascular
    "Hereditary Hemorrhagic Telangiectasia VCEP": "Vascular Disorders",
    # Epilepsy / Neurological
    "Epilepsy Sodium Channel VCEP": "Epilepsy / Neurological",
}

CATEGORY_ORDER = [
    "Cancer Predisposition", "Metabolic Disorders", "Cardiovascular Disorders", "Hematologic Disorders",
    "Immunodeficiency Disorders", "Neurodevelopmental Disorders", "Neuromuscular Disorders",
    "Ophthalmologic Disorders", "Hearing Disorders", "RASopathies",
    "Vascular Disorders", "Epilepsy / Neurological", "Other"
]


def get_category(expert_panel: str) -> str:
    """Map an Expert Panel (VCEP) to its broader disease category."""
    if not expert_panel or expert_panel == "N/A":
        return "Other"
    stripped = expert_panel.strip()
    category = VCEP_TO_CATEGORY.get(stripped)
    if category:
        return category
    # Case-insensitive fallback for minor formatting differences.
    lower = stripped.lower()
    for vcep, cat in VCEP_TO_CATEGORY.items():
        if vcep.lower() == lower:
            return cat
    return "Other"
