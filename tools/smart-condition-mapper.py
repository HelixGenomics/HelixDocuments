#!/usr/bin/env python3
"""Smart condition-to-PRS-trait mapper.

Handles free-text clinical conditions from PGP profiles and maps them to
our PRS scorecard trait names. Uses keyword extraction, abbreviation
expansion, and fuzzy matching to handle the wide variety of clinical naming.
"""

import re

# Our PRS traits (the ones we have PGS models for)
PRS_TRAITS = {
    # Cardiovascular
    'hypertension', 'coronary_artery_disease', 'myocardial_infarction',
    'atrial_fibrillation', 'heart_failure', 'stroke',
    'venous_thromboembolism', 'peripheral_artery_disease',
    # Metabolic
    'type_1_diabetes', 'type_2_diabetes', 'hypercholesterolemia',
    'hypertriglyceridemia', 'obesity', 'hypothyroidism', 'hyperthyroidism',
    'gout', 'metabolic_syndrome',
    # Cancer
    'breast_cancer', 'prostate_cancer', 'colorectal_cancer', 'lung_cancer',
    'melanoma', 'skin_cancer', 'basal_cell_carcinoma', 'thyroid_cancer',
    'bladder_cancer', 'ovarian_cancer', 'pancreatic_cancer',
    'endometrial_cancer', 'gastric_cancer', 'renal_cell_carcinoma',
    # Respiratory
    'asthma', 'copd', 'sleep_apnea',
    # Neurological/Psychiatric
    'depression', 'bipolar_disorder', 'schizophrenia', 'anxiety',
    'adhd', 'alzheimers_disease', 'parkinsons_disease', 'epilepsy',
    'migraine', 'autism_spectrum',
    # Autoimmune
    'rheumatoid_arthritis', 'systemic_lupus_erythematosus', 'celiac_disease',
    'crohns_disease', 'ulcerative_colitis', 'inflammatory_bowel_disease',
    'psoriasis', 'multiple_sclerosis', 'type_1_diabetes',
    'ankylosing_spondylitis',
    # GI
    'gerd', 'irritable_bowel_syndrome', 'gallstones', 'kidney_stones',
    # Eye
    'myopia', 'glaucoma', 'macular_degeneration', 'astigmatism',
    # Skin/Other
    'eczema', 'allergic_rhinitis', 'rosacea', 'seborrheic_dermatitis',
    'alopecia', 'vitiligo', 'osteoporosis', 'osteoarthritis',
    'chronic_kidney_disease', 'iron_deficiency_anemia',
    'benign_prostatic_hyperplasia',
}

# Qualifiers to strip from clinical text before matching
STRIP_QUALIFIERS = [
    r'\b(?:chronic|acute|severe|mild|moderate|bilateral|left|right|lt|rt)\b',
    r'\b(?:unspecified|other|primary|secondary|familial|hereditary|congenital)\b',
    r'\b(?:benign|malignant|recurrent|persistent|intermittent|progressive)\b',
    r'\b(?:early|late|adult|childhood|juvenile|senile|post|pre)\b',
    r'\b(?:upper|lower|posterior|anterior|medial|lateral)\b',
    r'\b(?:without|with|due to|associated|related|induced)\b',
    r'\b\d{4}\b',  # years like "2005"
    r'\b\d+/\d+\b',  # dates like "6/01"
    r'[,;].*$',  # strip everything after comma/semicolon (often secondary details)
]

# Abbreviation expansion
ABBREVIATIONS = {
    'cad': 'coronary artery disease',
    'mi': 'myocardial infarction',
    'afib': 'atrial fibrillation',
    'a-fib': 'atrial fibrillation',
    'chf': 'congestive heart failure',
    'dvt': 'deep vein thrombosis',
    'pe': 'pulmonary embolism',
    'vte': 'venous thromboembolism',
    'pad': 'peripheral artery disease',
    'pvd': 'peripheral vascular disease',
    'htn': 'hypertension',
    'dm': 'diabetes mellitus',
    't1d': 'type 1 diabetes',
    't2d': 'type 2 diabetes',
    'dm1': 'type 1 diabetes',
    'dm2': 'type 2 diabetes',
    'copd': 'chronic obstructive pulmonary disease',
    'osa': 'obstructive sleep apnea',
    'gerd': 'gastroesophageal reflux disease',
    'ibs': 'irritable bowel syndrome',
    'ibd': 'inflammatory bowel disease',
    'uc': 'ulcerative colitis',
    'ra': 'rheumatoid arthritis',
    'sle': 'systemic lupus erythematosus',
    'ms': 'multiple sclerosis',
    'as': 'ankylosing spondylitis',
    'adhd': 'attention deficit hyperactivity disorder',
    'asd': 'autism spectrum disorder',
    'mdd': 'major depressive disorder',
    'bph': 'benign prostatic hyperplasia',
    'ckd': 'chronic kidney disease',
    'uti': 'urinary tract infection',
    'bcc': 'basal cell carcinoma',
    'rcc': 'renal cell carcinoma',
    'ida': 'iron deficiency anemia',
}

# Keyword-to-trait mapping (if keyword appears in normalized text, maps to trait)
KEYWORD_MAP = {
    # Cardiovascular
    'hypertension': 'hypertension',
    'high blood pressure': 'hypertension',
    'elevated blood pressure': 'hypertension',
    'coronary artery disease': 'coronary_artery_disease',
    'coronary atherosclerosis': 'coronary_artery_disease',
    'ischemic heart': 'coronary_artery_disease',
    'angina': 'coronary_artery_disease',
    'myocardial infarction': 'myocardial_infarction',
    'heart attack': 'myocardial_infarction',
    'atrial fibrillation': 'atrial_fibrillation',
    'atrial flutter': 'atrial_fibrillation',  # close enough for PRS
    'heart failure': 'heart_failure',
    'cardiomyopathy': 'heart_failure',
    'stroke': 'stroke',
    'cerebrovascular': 'stroke',
    'transient ischemic': 'stroke',
    'deep vein thrombosis': 'venous_thromboembolism',
    'pulmonary embolism': 'venous_thromboembolism',
    'factor v leiden': 'venous_thromboembolism',
    'leiden factor': 'venous_thromboembolism',
    'thrombophilia': 'venous_thromboembolism',
    'peripheral vascular': 'peripheral_artery_disease',
    'peripheral artery': 'peripheral_artery_disease',
    'claudication': 'peripheral_artery_disease',

    # Metabolic
    'type 1 diabetes': 'type_1_diabetes',
    'type i diabetes': 'type_1_diabetes',
    'insulin dependent diabetes': 'type_1_diabetes',
    'type 2 diabetes': 'type_2_diabetes',
    'type ii diabetes': 'type_2_diabetes',
    'non-insulin dependent diabetes': 'type_2_diabetes',
    'diabetes mellitus': 'type_2_diabetes',  # default to T2D
    'diabetes type 2': 'type_2_diabetes',
    'diabetes, type 2': 'type_2_diabetes',
    'diabetes type 1': 'type_1_diabetes',
    'diabetes, type 1': 'type_1_diabetes',
    'hypercholesterolemia': 'hypercholesterolemia',
    'high cholesterol': 'hypercholesterolemia',
    'elevated cholesterol': 'hypercholesterolemia',
    'hyperlipidemia': 'hypercholesterolemia',
    'dyslipidemia': 'hypercholesterolemia',
    'hypertriglyceridemia': 'hypertriglyceridemia',
    'high triglycerides': 'hypertriglyceridemia',
    'elevated triglycerides': 'hypertriglyceridemia',
    'obesity': 'obesity',
    'morbid obesity': 'obesity',
    'hypothyroidism': 'hypothyroidism',
    'hashimoto': 'hypothyroidism',
    'underactive thyroid': 'hypothyroidism',
    'hyperthyroidism': 'hyperthyroidism',
    'graves disease': 'hyperthyroidism',
    'overactive thyroid': 'hyperthyroidism',
    'gout': 'gout',
    'hyperuricemia': 'gout',
    'homocystinemia': 'hypercholesterolemia',  # related metabolic risk

    # Cancer
    'breast cancer': 'breast_cancer',
    'breast carcinoma': 'breast_cancer',
    'prostate cancer': 'prostate_cancer',
    'prostate carcinoma': 'prostate_cancer',
    'colorectal cancer': 'colorectal_cancer',
    'colon cancer': 'colorectal_cancer',
    'rectal cancer': 'colorectal_cancer',
    'bowel cancer': 'colorectal_cancer',
    'lung cancer': 'lung_cancer',
    'pulmonary neoplasm': 'lung_cancer',
    'melanoma': 'melanoma',
    'skin cancer': 'skin_cancer',
    'basal cell carcinoma': 'basal_cell_carcinoma',
    'squamous cell carcinoma': 'skin_cancer',
    'thyroid cancer': 'thyroid_cancer',
    'thyroid carcinoma': 'thyroid_cancer',
    'bladder cancer': 'bladder_cancer',
    'ovarian cancer': 'ovarian_cancer',
    'pancreatic cancer': 'pancreatic_cancer',
    'endometrial cancer': 'endometrial_cancer',
    'uterine cancer': 'endometrial_cancer',
    'gastric cancer': 'gastric_cancer',
    'stomach cancer': 'gastric_cancer',
    'renal cell carcinoma': 'renal_cell_carcinoma',
    'kidney cancer': 'renal_cell_carcinoma',

    # Respiratory
    'asthma': 'asthma',
    'reactive airway': 'asthma',
    'obstructive pulmonary': 'copd',
    'emphysema': 'copd',
    'chronic bronchitis': 'copd',
    'sleep apnea': 'sleep_apnea',

    # Neurological/Psychiatric
    'depression': 'depression',
    'depressive disorder': 'depression',
    'bipolar': 'bipolar_disorder',
    'manic depression': 'bipolar_disorder',
    'schizophrenia': 'schizophrenia',
    'anxiety': 'anxiety',
    'panic disorder': 'anxiety',
    'generalized anxiety': 'anxiety',
    'attention deficit': 'adhd',
    'alzheimer': 'alzheimers_disease',
    'dementia': 'alzheimers_disease',
    'parkinson': 'parkinsons_disease',
    'epilepsy': 'epilepsy',
    'seizure disorder': 'epilepsy',
    'migraine': 'migraine',
    'autism': 'autism_spectrum',
    'asperger': 'autism_spectrum',

    # Autoimmune
    'rheumatoid arthritis': 'rheumatoid_arthritis',
    'lupus': 'systemic_lupus_erythematosus',
    'celiac': 'celiac_disease',
    'crohn': 'crohns_disease',
    'ulcerative colitis': 'ulcerative_colitis',
    'psoriasis': 'psoriasis',
    'psoriatic': 'psoriasis',
    'multiple sclerosis': 'multiple_sclerosis',
    'ankylosing spondylitis': 'ankylosing_spondylitis',

    # GI
    'gastroesophageal reflux': 'gerd',
    'acid reflux': 'gerd',
    'irritable bowel': 'irritable_bowel_syndrome',
    'cholelithiasis': 'gallstones',
    'gallstones': 'gallstones',
    'gallbladder stones': 'gallstones',
    'kidney stones': 'kidney_stones',
    'nephrolithiasis': 'kidney_stones',
    'renal calculi': 'kidney_stones',

    # Eye
    'myopia': 'myopia',
    'nearsightedness': 'myopia',
    'near-sightedness': 'myopia',
    'glaucoma': 'glaucoma',
    'macular degeneration': 'macular_degeneration',
    'astigmatism': 'astigmatism',
    'farsightedness': None,  # no PGS model typically

    # Skin/Other
    'eczema': 'eczema',
    'atopic dermatitis': 'eczema',
    'allergic rhinitis': 'allergic_rhinitis',
    'hay fever': 'allergic_rhinitis',
    'rosacea': 'rosacea',
    'seborrheic dermatitis': 'seborrheic_dermatitis',
    'dandruff': 'seborrheic_dermatitis',
    'alopecia': 'alopecia',
    'pattern baldness': 'alopecia',
    'hair loss': 'alopecia',
    'vitiligo': 'vitiligo',
    'osteoporosis': 'osteoporosis',
    'osteoarthritis': 'osteoarthritis',
    'degenerative joint': 'osteoarthritis',
    'chronic kidney disease': 'chronic_kidney_disease',
    'renal insufficiency': 'chronic_kidney_disease',
    'iron deficiency anemia': 'iron_deficiency_anemia',
    'prostatic hyperplasia': 'benign_prostatic_hyperplasia',
    'enlarged prostate': 'benign_prostatic_hyperplasia',
}


def normalize_condition(text):
    """Normalize a clinical condition text for matching."""
    # Lowercase
    text = text.lower().strip()

    # Expand known abbreviations (check whole text for abbreviation-only entries)
    words = text.split()
    if len(words) <= 3:
        for abbr, expansion in ABBREVIATIONS.items():
            if text.strip().replace('.', '') == abbr:
                text = expansion
                break

    # Strip qualifiers
    for pattern in STRIP_QUALIFIERS:
        text = re.sub(pattern, ' ', text, flags=re.IGNORECASE)

    # Normalize whitespace
    text = re.sub(r'\s+', ' ', text).strip()

    return text


def map_condition_to_trait(raw_condition):
    """Map a clinical condition to a PRS trait name.

    Returns (trait_name, confidence) where:
    - trait_name is the standardized PRS trait or None
    - confidence is 'high' (exact keyword match), 'medium' (after normalization), or 'low' (partial)
    """
    raw_lower = raw_condition.lower().strip()

    # 1. Direct keyword match on raw text (highest confidence)
    for keyword, trait in KEYWORD_MAP.items():
        if keyword in raw_lower:
            if trait is not None:
                return trait, 'high'
            else:
                return None, 'skip'  # Known non-PRS condition

    # 2. Abbreviation expansion — check each word
    words = raw_lower.replace('.', '').replace(',', ' ').split()
    expanded_words = []
    for w in words:
        w_clean = w.strip()
        if w_clean in ABBREVIATIONS:
            expanded_words.append(ABBREVIATIONS[w_clean])
        else:
            expanded_words.append(w_clean)
    expanded = ' '.join(expanded_words)
    if expanded != raw_lower:
        for keyword, trait in KEYWORD_MAP.items():
            if keyword in expanded:
                if trait is not None:
                    return trait, 'high'

    # 3. Normalize and try again
    normalized = normalize_condition(raw_condition)
    if normalized != raw_lower:
        for keyword, trait in KEYWORD_MAP.items():
            if keyword in normalized:
                if trait is not None:
                    return trait, 'medium'

    # 4. Partial word matching (lower confidence)
    # Check if any key PRS trait word appears
    for trait in PRS_TRAITS:
        trait_words = trait.replace('_', ' ').split()
        if len(trait_words) == 1:
            if trait_words[0] in raw_lower.split():
                return trait, 'low'
        elif all(w in raw_lower for w in trait_words):
            return trait, 'medium'

    return None, 'none'


def test_mapper():
    """Test with the example conditions from PGP profiles."""
    test_conditions = [
        "Chicken Pox",
        "Measles",
        "Vental Hernia",
        "Elective Breast Implants, Bilateral",
        "Full Tear Left Rotator Cuff",
        "menorrhagia",
        "Left wrist fracture",
        "Cholelithiasis without obstruction",
        "Environmental allergies",
        "Chronic Sinusitis; Left deviated septum",
        "Pneumonia and influenza",
        "Post herpetic Syndrome - left lower face",
        "Osteoarthritis, Bilateral Knees 2005",
        "Positive Leiden factor V (heterozygous presence only) 6/01",
        "Atrial Flutter",
        "Rt CAD",
        "Familial hypertriglyceridemia",
        "Diabetes, Type 2",
        "osteoarthritis bilateral hips",
        "Homocystinemia",
        "Chest Pain, Acute hypertension",
        "Nearsightedness",
        "Chronic hypertension",
        "Glaucoma",
        "Farsightedness",
        "Dry macular degeneration",
        "ASTHMA",
        "OTHER UNSPECIFIED HYPERLIPIDEMIA",
        "UNSPECIFIED SLEEP APNEA",
        "SEBORRHEIC DERMATITIS, UNSPECIFIED",
        "High cholesterol (hypercholesterolemia)",
        "Myopia (Nearsightedness)",
        "Type 2 Diabetes",
        "High blood pressure",
        "Irritable bowel syndrome (IBS)",
    ]

    print(f"{'Condition':<55} {'Trait':<30} {'Confidence'}")
    print('=' * 100)
    mapped = 0
    for cond in test_conditions:
        trait, conf = map_condition_to_trait(cond)
        if trait:
            mapped += 1
        marker = '*' if conf in ('high', 'medium') else (' ' if conf == 'low' else '.')
        trait_str = trait or '-'
        print(f"{marker} {cond:<53} {trait_str:<30} {conf}")

    print(f'\n{mapped}/{len(test_conditions)} conditions mapped to PRS traits')


if __name__ == '__main__':
    test_mapper()
