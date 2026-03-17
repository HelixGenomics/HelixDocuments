#!/usr/bin/env python3
"""Deep scrape PGP profiles for survey responses + medical records.

Re-scrapes profiles that have genotype data, extracting the full survey
Q&A pairs which contain structured disease diagnoses, demographics, and
self-reported traits. This is the real phenotype ground truth for PRS validation.
"""

import os, json, time, re
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from datetime import datetime

OUTPUT_DIR = '/workspace/pgp-pipeline'
PHENOTYPES_FILE = f'{OUTPUT_DIR}/pgp-phenotypes.json'
OUTPUT = f'{OUTPUT_DIR}/pgp-deep-phenotypes.json'

# Disease survey categories and their typical question patterns
DISEASE_CATEGORIES = [
    'Cancers', 'Endocrine', 'Metabolic', 'Blood', 'Nervous System',
    'Vision', 'Hearing', 'Circulatory', 'Respiratory', 'Digestive',
    'Genitourinary', 'Skin', 'Musculoskeletal', 'Congenital',
]

# Map of common PGP survey conditions to PGS-relevant traits
# This lets us match self-reported diagnoses to our PRS models
CONDITION_TO_TRAIT = {
    # Cardiovascular
    'hypertension': 'hypertension',
    'high blood pressure': 'hypertension',
    'coronary artery disease': 'coronary_artery_disease',
    'heart attack': 'myocardial_infarction',
    'myocardial infarction': 'myocardial_infarction',
    'atrial fibrillation': 'atrial_fibrillation',
    'heart failure': 'heart_failure',
    'stroke': 'stroke',
    'deep vein thrombosis': 'venous_thromboembolism',
    'pulmonary embolism': 'venous_thromboembolism',

    # Metabolic
    'type 1 diabetes': 'type_1_diabetes',
    'type 2 diabetes': 'type_2_diabetes',
    'diabetes': 'type_2_diabetes',
    'high cholesterol': 'hypercholesterolemia',
    'hypercholesterolemia': 'hypercholesterolemia',
    'high triglycerides': 'hypertriglyceridemia',
    'hypertriglyceridemia': 'hypertriglyceridemia',
    'obesity': 'obesity',
    'hypothyroidism': 'hypothyroidism',
    'hyperthyroidism': 'hyperthyroidism',
    'gout': 'gout',

    # Cancer
    'breast cancer': 'breast_cancer',
    'prostate cancer': 'prostate_cancer',
    'colorectal cancer': 'colorectal_cancer',
    'colon cancer': 'colorectal_cancer',
    'lung cancer': 'lung_cancer',
    'melanoma': 'melanoma',
    'skin cancer': 'skin_cancer',
    'basal cell carcinoma': 'basal_cell_carcinoma',
    'thyroid cancer': 'thyroid_cancer',
    'bladder cancer': 'bladder_cancer',
    'ovarian cancer': 'ovarian_cancer',

    # Respiratory
    'asthma': 'asthma',
    'asthma (adult)': 'asthma',
    'asthma (childhood)': 'asthma',
    'copd': 'copd',
    'chronic obstructive pulmonary disease': 'copd',
    'sleep apnea': 'sleep_apnea',

    # Neurological / Psychiatric
    'depression': 'depression',
    'major depression': 'depression',
    'bipolar disorder': 'bipolar_disorder',
    'schizophrenia': 'schizophrenia',
    'anxiety': 'anxiety',
    'adhd': 'adhd',
    'attention deficit hyperactivity disorder': 'adhd',
    'alzheimer': 'alzheimers_disease',
    'parkinson': 'parkinsons_disease',
    'epilepsy': 'epilepsy',
    'migraine': 'migraine',

    # Autoimmune
    'rheumatoid arthritis': 'rheumatoid_arthritis',
    'lupus': 'systemic_lupus_erythematosus',
    'celiac disease': 'celiac_disease',
    'crohn': 'crohns_disease',
    'ulcerative colitis': 'ulcerative_colitis',
    'inflammatory bowel disease': 'inflammatory_bowel_disease',
    'psoriasis': 'psoriasis',
    'multiple sclerosis': 'multiple_sclerosis',

    # GI
    'gastroesophageal reflux': 'gerd',
    'gerd': 'gerd',
    'acid reflux': 'gerd',
    'irritable bowel syndrome': 'irritable_bowel_syndrome',
    'ibs': 'irritable_bowel_syndrome',
    'gallstones': 'gallstones',
    'kidney stones': 'kidney_stones',

    # Eye
    'myopia': 'myopia',
    'nearsightedness': 'myopia',
    'glaucoma': 'glaucoma',
    'macular degeneration': 'macular_degeneration',
    'cystoid macular degeneration': 'macular_degeneration',
    'astigmatism': 'astigmatism',
    'regular astigmatism': 'astigmatism',

    # Other
    'osteoporosis': 'osteoporosis',
    'osteoarthritis': 'osteoarthritis',
    'eczema': 'eczema',
    'allergic rhinitis': 'allergic_rhinitis',
    'hay fever': 'allergic_rhinitis',
    'rosacea': 'rosacea',
    'seborrheic dermatitis': 'seborrheic_dermatitis',
    'alopecia': 'alopecia',
    'vitiligo': 'vitiligo',

    # ICD-style clinical condition names (from PGP medical records)
    'unspecified essential hypertension': 'hypertension',
    'essential hypertension': 'hypertension',
    'hyperlipidemia': 'hypercholesterolemia',
    'unspecified hyperlipidemia': 'hypercholesterolemia',
    'other unspecified hyperlipidemia': 'hypercholesterolemia',
    'pure hypercholesterolemia': 'hypercholesterolemia',
    'unspecified sleep apnea': 'sleep_apnea',
    'obstructive sleep apnea': 'sleep_apnea',
    'diabetes mellitus': 'type_2_diabetes',
    'non-insulin dependent diabetes': 'type_2_diabetes',
    'insulin dependent diabetes': 'type_1_diabetes',
    'depressive disorder': 'depression',
    'major depressive disorder': 'depression',
    'generalized anxiety disorder': 'anxiety',
    'anxiety disorder': 'anxiety',
    'panic disorder': 'anxiety',
    'malignant neoplasm of breast': 'breast_cancer',
    'malignant neoplasm of prostate': 'prostate_cancer',
    'malignant neoplasm of colon': 'colorectal_cancer',
    'malignant melanoma': 'melanoma',
    'atopic dermatitis': 'eczema',
    'allergic rhinitis due to pollen': 'allergic_rhinitis',
    'chronic kidney disease': 'chronic_kidney_disease',
    'iron deficiency anemia': 'iron_deficiency_anemia',
    'hyperuricemia': 'gout',
    'benign prostatic hyperplasia': 'benign_prostatic_hyperplasia',
    'chronic obstructive airway disease': 'copd',
    'unspecified asthma': 'asthma',
    'extrinsic asthma': 'asthma',
    'intrinsic asthma': 'asthma',
    'coronary atherosclerosis': 'coronary_artery_disease',
    'cerebrovascular disease': 'stroke',
    'peripheral vascular disease': 'peripheral_artery_disease',
}

def ts():
    return datetime.now().strftime('%H:%M:%S')


def extract_surveys(html):
    """Extract all survey Q&A pairs from profile HTML."""
    # Survey responses are in: <tr class="survey_result_NNN ui-helper-hidden"><td>Q</td><td class="hoverable">A</td></tr>
    pattern = r'<tr class="survey_result_(\d+) ui-helper-hidden">\s*<td>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*</tr>'
    matches = re.findall(pattern, html, re.DOTALL)

    # Find survey names mapped to IDs
    name_pattern = r'(PGP [^<]+?(?:Survey|Assessment)[^<]*?)\s*</th>.*?survey_result_(\d+)'
    name_map = {}
    for m in re.finditer(name_pattern, html, re.DOTALL):
        name = m.group(1).replace('&amp;', '&').strip()
        sid = m.group(2)
        if sid not in name_map:
            name_map[sid] = name

    # Group by survey ID
    surveys = {}
    for sid, q, a in matches:
        q = re.sub(r'<[^>]+>', '', q).strip()
        a = re.sub(r'<[^>]+>', '', a).strip()
        if q and a:
            if sid not in surveys:
                surveys[sid] = {'name': name_map.get(sid, f'Survey #{sid}'), 'responses': []}
            surveys[sid]['responses'].append({'question': q, 'answer': a})

    return surveys


def extract_diagnoses(surveys):
    """Extract structured diagnoses from survey responses."""
    diagnoses = []
    seen = set()

    for sid, survey in surveys.items():
        survey_name = survey['name']
        for resp in survey['responses']:
            q = resp['question'].lower()
            a = resp['answer']

            # Pattern 1: "Have you ever been diagnosed with any of the following conditions?" -> comma list
            if 'diagnosed' in q and 'following' in q and a and a.lower() not in ('', 'no', 'none', 'no response', 'timestamp'):
                # Split comma-separated conditions
                conditions = [c.strip() for c in a.split(',')]
                for cond in conditions:
                    cond_clean = cond.strip()
                    if cond_clean and cond_clean.lower() not in ('no', 'none'):
                        key = cond_clean.lower()
                        if key not in seen:
                            seen.add(key)
                            # Try to map to standard trait
                            trait = None
                            for pattern, t in CONDITION_TO_TRAIT.items():
                                if pattern in key:
                                    trait = t
                                    break
                            diagnoses.append({
                                'condition': cond_clean,
                                'trait': trait,
                                'source': survey_name,
                                'type': 'self_reported_list',
                            })

            # Pattern 2: "Have you ever been diagnosed with...? [Condition Name]" -> Yes/No
            elif 'diagnosed' in q and '[' in resp['question']:
                bracket_match = re.search(r'\[([^\]]+)\]', resp['question'])
                if bracket_match and a.lower() in ('yes', 'true'):
                    cond_clean = bracket_match.group(1).strip()
                    key = cond_clean.lower()
                    if key not in seen:
                        seen.add(key)
                        trait = None
                        for pattern, t in CONDITION_TO_TRAIT.items():
                            if pattern in key:
                                trait = t
                                break
                        diagnoses.append({
                            'condition': cond_clean,
                            'trait': trait,
                            'source': survey_name,
                            'type': 'self_reported_yn',
                        })

    return diagnoses


def extract_demographics(surveys):
    """Extract demographic info from survey responses."""
    demo = {}

    for sid, survey in surveys.items():
        for resp in survey['responses']:
            q = resp['question'].lower()
            a = resp['answer']

            if not a or a.lower() in ('no response', ''):
                continue

            if 'year of birth' in q and not demo.get('birth_year'):
                m = re.search(r'(\d{4})', a)
                if m:
                    demo['birth_year'] = int(m.group(1))
                elif 'years' in a.lower():
                    # "30-39 years" format
                    demo['birth_year_range'] = a

            if 'sex' in q and 'gender' in q and not demo.get('sex'):
                demo['sex'] = a

            if ('anatomical sex' in q or 'sex at birth' in q) and not demo.get('sex_at_birth'):
                demo['sex_at_birth'] = a

            if 'race' in q and 'ethnicity' not in q and not demo.get('race'):
                demo['race'] = a

            if 'ethnicity' in q and 'race' not in q and not demo.get('ethnicity'):
                demo['ethnicity'] = a

            if 'race' in q and 'ethnicity' in q and not demo.get('race_ethnicity'):
                demo['race_ethnicity'] = a

            if 'grandmother' in q and 'country' in q:
                if 'maternal' in q:
                    demo['maternal_grandmother_origin'] = a
                elif 'paternal' in q:
                    demo['paternal_grandmother_origin'] = a

            if 'grandfather' in q and 'country' in q:
                if 'maternal' in q:
                    demo['maternal_grandfather_origin'] = a
                elif 'paternal' in q:
                    demo['paternal_grandfather_origin'] = a

            if 'month of birth' in q and not demo.get('birth_month'):
                demo['birth_month'] = a

            if 'what is your age' in q and not demo.get('reported_age'):
                m = re.search(r'(\d+)', a)
                if m:
                    demo['reported_age'] = int(m.group(1))

    return demo


def extract_profile_basics(html):
    """Extract basic profile data (height, weight, DOB, conditions, tests)."""
    data = {
        'height_cm': None, 'weight_kg': None, 'dob': None, 'age': None,
        'gender': None, 'race': None, 'blood_type': None,
        'conditions': [], 'test_results': [], 'procedures': [],
    }

    # Height
    m = re.search(r'(\d+)\s*(?:ft|\')\s*(\d+)\s*(?:in|\")\s*\((\d+)\s*cm\)', html)
    if m:
        data['height_cm'] = int(m.group(3))

    # Weight
    m = re.search(r'(\d+)\s*lbs?\s*\((\d+)\s*kg\)', html)
    if m:
        data['weight_kg'] = float(m.group(2))

    # DOB
    m = re.search(r'(\d{4}-\d{2}-\d{2})\s*\((\d+) years old\)', html)
    if m:
        data['dob'] = m.group(1)
        data['age'] = int(m.group(2))

    # Gender
    m = re.search(r'Gender\s*</td>\s*<td[^>]*>\s*([^<]+)', html)
    if m and m.group(1).strip():
        data['gender'] = m.group(1).strip()

    # Race
    m = re.search(r'Race(?:/ethnicity)?\s*</td>\s*<td[^>]*>\s*([^<]+)', html)
    if m and m.group(1).strip():
        data['race'] = m.group(1).strip()

    # Blood type
    m = re.search(r'Blood Type\s*</td>\s*<td[^>]*>\s*([^<]+)', html)
    if m and m.group(1).strip() in ('A+', 'A-', 'B+', 'B-', 'AB+', 'AB-', 'O+', 'O-'):
        data['blood_type'] = m.group(1).strip()

    # Conditions table
    cond_block = re.search(r'Conditions</h\d>.*?<table[^>]*>(.*?)</table>', html, re.DOTALL)
    if cond_block:
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', cond_block.group(1), re.DOTALL)
        for row in rows:
            cells = re.findall(r'<td[^>]*>(.*?)</td>', row, re.DOTALL)
            if cells:
                name = re.sub(r'<[^>]+>', '', cells[0]).strip()
                if name and name not in ('Name', 'Start Date'):
                    data['conditions'].append(name)

    # Test results table
    test_block = re.search(r'Test Results</h\d>.*?<table[^>]*>(.*?)</table>', html, re.DOTALL)
    if test_block:
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', test_block.group(1), re.DOTALL)
        for row in rows:
            cells = re.findall(r'<td[^>]*>(.*?)</td>', row, re.DOTALL)
            if len(cells) >= 2:
                name = re.sub(r'<[^>]+>', '', cells[0]).strip()
                result = re.sub(r'<[^>]+>', '', cells[1]).strip()
                if name and name not in ('Name', 'Result'):
                    data['test_results'].append({'name': name, 'result': result})

    # Procedures table
    proc_block = re.search(r'Procedures</h\d>.*?<table[^>]*>(.*?)</table>', html, re.DOTALL)
    if proc_block:
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', proc_block.group(1), re.DOTALL)
        for row in rows:
            cells = re.findall(r'<td[^>]*>(.*?)</td>', row, re.DOTALL)
            if cells:
                name = re.sub(r'<[^>]+>', '', cells[0]).strip()
                if name and name not in ('Name', 'Date'):
                    data['procedures'].append(name)

    # 23andMe file links
    file_matches = re.findall(r'/user_file/download/(\d+)', html)
    data['genotype_files'] = list(set(file_matches))

    return data


def scrape_one_deep(hex_id):
    """Deep scrape a single profile."""
    url = f'https://my.pgp-hms.org/profile_public?hex={hex_id}'

    for attempt in range(3):
        try:
            req = Request(url, headers={'User-Agent': 'HelixGenomics/1.0 (research)'})
            resp = urlopen(req, timeout=30)
            html = resp.read().decode('utf-8', errors='replace')

            # Basic profile data
            basics = extract_profile_basics(html)

            # Full survey extraction
            surveys = extract_surveys(html)

            # Structured diagnoses from surveys
            diagnoses = extract_diagnoses(surveys)

            # Demographics from surveys
            survey_demographics = extract_demographics(surveys)

            # Map clinical conditions from profile to traits
            clinical_diagnoses = []
            for c in basics.get('conditions', []):
                cond_lower = c.lower().strip()
                trait = None
                for pattern, t in CONDITION_TO_TRAIT.items():
                    if pattern in cond_lower:
                        trait = t
                        break
                clinical_diagnoses.append({
                    'condition': c,
                    'trait': trait,
                    'source': 'medical_records',
                    'type': 'clinical_diagnosis',
                })

            # Merge survey diagnoses with clinical diagnoses
            all_diagnoses = clinical_diagnoses + diagnoses
            all_conditions = set()
            for c in basics.get('conditions', []):
                all_conditions.add(c.lower().strip())
            for d in diagnoses:
                all_conditions.add(d['condition'].lower().strip())

            # Count survey responses
            total_qa = sum(len(s['responses']) for s in surveys.values())
            n_surveys = len(surveys)

            # Build result
            result = {
                'hex_id': hex_id,
                'profile_url': url,
                # Basic measurements
                'height_cm': basics['height_cm'],
                'weight_kg': basics['weight_kg'],
                'dob': basics['dob'],
                'age': basics['age'],
                'gender': basics['gender'],
                'race': basics['race'],
                'blood_type': basics['blood_type'],
                # Survey demographics (often richer)
                'survey_demographics': survey_demographics,
                # Medical records
                'clinical_conditions': basics['conditions'],  # From medical records table (ICD-coded)
                'clinical_diagnoses': clinical_diagnoses,  # Mapped to traits
                'survey_diagnoses': diagnoses,  # From survey self-reports
                'all_diagnoses': all_diagnoses,  # Combined
                'all_conditions': sorted(all_conditions),
                'n_conditions': len(all_conditions),
                'test_results': basics['test_results'],
                'procedures': basics['procedures'],
                # Survey metadata
                'n_surveys': n_surveys,
                'n_survey_qa': total_qa,
                'survey_names': [s['name'] for s in surveys.values()],
                # Raw survey data (for detailed analysis)
                'surveys_raw': {sid: s for sid, s in surveys.items()},
                # File matching
                'genotype_files': basics['genotype_files'],
            }

            # Compute rich validation score
            score = 0
            if basics['height_cm']: score += 5
            if basics['weight_kg']: score += 5
            if basics['gender'] or survey_demographics.get('sex_at_birth'): score += 2
            if basics['dob'] or survey_demographics.get('birth_year'): score += 2
            if basics['race'] or survey_demographics.get('race'): score += 2
            if basics['blood_type']: score += 2
            # Clinical diagnoses are most valuable (doctor-confirmed)
            score += len([d for d in clinical_diagnoses if d.get('trait')]) * 5
            score += len(clinical_diagnoses) * 2
            # Survey self-reports are useful but less reliable
            score += len([d for d in diagnoses if d.get('trait')]) * 3
            score += len(all_conditions) * 1
            score += len(basics['test_results'])
            score += n_surveys * 1
            score += min(total_qa, 50)  # Cap survey Q&A contribution
            result['validation_score'] = score

            # Traits we can validate PRS against (from both clinical + survey)
            validatable_traits = list(set(d['trait'] for d in all_diagnoses if d.get('trait')))
            result['validatable_traits'] = sorted(validatable_traits)
            result['n_validatable_traits'] = len(validatable_traits)

            return result

        except HTTPError as e:
            if e.code == 429:
                time.sleep(5)
                continue
            elif e.code == 404:
                return {'hex_id': hex_id, 'error': 'not_found', 'validation_score': 0}
            time.sleep(2)
        except Exception as e:
            time.sleep(2)

    return {'hex_id': hex_id, 'error': 'failed', 'validation_score': 0}


def main():
    print(f'[{ts()}] Deep PGP Phenotype Scraper', flush=True)

    # Load initial scrape results to get matched profiles
    if os.path.exists(PHENOTYPES_FILE):
        with open(PHENOTYPES_FILE) as f:
            initial = json.load(f)
        # Get hex IDs that have genotype files (matched to our downloads)
        matched = [p['hex_id'] for p in initial['profiles']
                   if not p.get('error') and p.get('matched_files')]
        # Also include high-richness profiles even if not matched
        rich = [p['hex_id'] for p in initial['profiles']
                if not p.get('error') and p.get('richness_score', 0) > 5
                and p['hex_id'] not in set(matched)]
        hex_ids = matched + rich
        print(f'  {len(matched)} matched to genotypes + {len(rich)} rich profiles = {len(hex_ids)} total', flush=True)
    else:
        # Fallback: scrape from genetic data page
        print(f'  No initial phenotypes file, scraping from genetic data page...', flush=True)
        url = 'https://my.pgp-hms.org/public_genetic_data?data_type=23andMe'
        req = Request(url, headers={'User-Agent': 'HelixGenomics/1.0 (research)'})
        resp = urlopen(req, timeout=60)
        html = resp.read().decode('utf-8', errors='replace')
        hex_ids = list(set(re.findall(r'profile/([a-zA-Z0-9]+)', html)))
        print(f'  Found {len(hex_ids)} participants', flush=True)

    # Deep scrape
    print(f'\n[{ts()}] Deep scraping {len(hex_ids)} profiles...', flush=True)

    profiles = []
    done = 0
    errors = 0

    with ThreadPoolExecutor(max_workers=10) as pool:
        futures = {pool.submit(scrape_one_deep, h): h for h in hex_ids}
        for future in as_completed(futures):
            result = future.result()
            profiles.append(result)
            done += 1
            if result.get('error'):
                errors += 1

            if done % 50 == 0:
                print(f'  [{ts()}] Deep scraped {done}/{len(hex_ids)} ({errors} errors)', flush=True)

            time.sleep(0.15)

    # Sort by validation score
    profiles.sort(key=lambda p: p.get('validation_score', 0), reverse=True)

    # Analysis
    valid = [p for p in profiles if not p.get('error')]
    with_height = [p for p in valid if p.get('height_cm')]
    with_weight = [p for p in valid if p.get('weight_kg')]
    with_clinical = [p for p in valid if p.get('clinical_conditions')]
    with_survey_diag = [p for p in valid if p.get('survey_diagnoses')]
    with_validatable = [p for p in valid if p.get('n_validatable_traits', 0) > 0]
    with_surveys = [p for p in valid if p.get('n_surveys', 0) > 0]

    # Trait coverage for PRS validation
    trait_counts = {}
    for p in valid:
        for t in p.get('validatable_traits', []):
            trait_counts[t] = trait_counts.get(t, 0) + 1

    # All conditions frequency
    all_cond_counts = {}
    for p in valid:
        for c in p.get('all_conditions', []):
            all_cond_counts[c] = all_cond_counts.get(c, 0) + 1

    print(f'\n{"="*60}', flush=True)
    print(f'Deep Phenotype Scraping Results', flush=True)
    print(f'{"="*60}', flush=True)
    print(f'Total scraped: {len(profiles)}', flush=True)
    print(f'Valid: {len(valid)}, Errors: {errors}', flush=True)
    print(f'', flush=True)
    print(f'Phenotype Coverage:', flush=True)
    print(f'  Height: {len(with_height)}', flush=True)
    print(f'  Weight: {len(with_weight)}', flush=True)
    print(f'  Survey responses: {len(with_surveys)}', flush=True)
    print(f'  Clinical diagnoses (medical records): {len(with_clinical)}', flush=True)
    print(f'  Survey diagnoses (self-reported): {len(with_survey_diag)}', flush=True)
    print(f'  Has PRS-validatable traits: {len(with_validatable)}', flush=True)

    print(f'\nPRS-Validatable Traits (self-reported diagnoses matched to PGS models):', flush=True)
    for trait, count in sorted(trait_counts.items(), key=lambda x: -x[1]):
        print(f'  {trait}: {count} participants', flush=True)

    print(f'\nTop 30 All Conditions (including unmapped):', flush=True)
    for cond, count in sorted(all_cond_counts.items(), key=lambda x: -x[1])[:30]:
        print(f'  {cond}: {count}', flush=True)

    print(f'\nTop 15 Best Profiles for PRS Validation:', flush=True)
    for p in profiles[:15]:
        hid = p.get('hex_id', '?')
        vs = p.get('validation_score', 0)
        h = f"{p.get('height_cm','')}cm" if p.get('height_cm') else '-'
        w = f"{p.get('weight_kg','')}kg" if p.get('weight_kg') else '-'
        nc = p.get('n_conditions', 0)
        nv = p.get('n_validatable_traits', 0)
        ns = p.get('n_surveys', 0)
        nqa = p.get('n_survey_qa', 0)
        print(f'  {hid}: score={vs}, height={h}, weight={w}, {nc} conditions, {nv} validatable traits, {ns} surveys ({nqa} Q&A)', flush=True)

    # Save results - strip raw survey data for profiles with low scores to save space
    for p in profiles:
        if p.get('validation_score', 0) < 10:
            p.pop('surveys_raw', None)

    with open(OUTPUT, 'w') as f:
        json.dump({
            'scraped_at': datetime.now().isoformat(),
            'total_profiles': len(profiles),
            'valid_profiles': len(valid),
            'summary': {
                'with_height': len(with_height),
                'with_weight': len(with_weight),
                'with_surveys': len(with_surveys),
                'with_clinical_diagnoses': len(with_clinical),
                'with_survey_diagnoses': len(with_survey_diag),
                'with_validatable_traits': len(with_validatable),
                'trait_counts': trait_counts,
                'top_conditions': dict(sorted(all_cond_counts.items(), key=lambda x: -x[1])[:100]),
            },
            'profiles': profiles,
        }, f, indent=2)

    print(f'\nSaved to {OUTPUT}', flush=True)
    print(f'[{ts()}] Done.', flush=True)


if __name__ == '__main__':
    main()
