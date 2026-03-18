#!/usr/bin/env python3
"""
Generate plain-English descriptions for PGS traits using Claude CLI.
Batches ~85 traits per call, saves results to pgs-trait-descriptions.json.
"""
import json, subprocess, sys, os, time, re

BATCH_SIZE = 40
OUTPUT_FILE = "/opt/helix/data/pgs-trait-descriptions.json"
TRAITS_FILE = "/tmp/unique-pgs-traits.json"

with open(TRAITS_FILE) as f:
    all_traits = json.load(f)

# Load existing descriptions if any (resume support)
existing = {}
if os.path.exists(OUTPUT_FILE):
    try:
        with open(OUTPUT_FILE) as f:
            existing = json.load(f)
        print(f"Loaded {len(existing)} existing descriptions")
    except:
        pass

# Filter out already-described traits
remaining = [t for t in all_traits if t not in existing]
print(f"Total traits: {len(all_traits)}, remaining: {len(remaining)}")

if not remaining:
    print("All traits already described!")
    sys.exit(0)

batches = [remaining[i:i+BATCH_SIZE] for i in range(0, len(remaining), BATCH_SIZE)]
print(f"Will process {len(batches)} batches of ~{BATCH_SIZE} traits each")

for batch_idx, batch in enumerate(batches):
    print(f"\n=== Batch {batch_idx+1}/{len(batches)} ({len(batch)} traits) ===")

    trait_list = "\n".join(f"- {t}" for t in batch)

    prompt = f"""You are writing short descriptions of genetic traits for a consumer DNA health report.
The audience is normal people with NO medical or genetics background.

For each trait below, write a 1-2 sentence plain-English description that explains:
1. What this trait actually means in simple words
2. What it means if someone scores high or low

Rules:
- Use simple language a 16 year old would understand
- No jargon - if you must use a medical term, explain it in parentheses
- Keep each description under 40 words
- Be direct and practical
- For obscure biomarkers, explain what they measure and why it matters
- For diseases, explain what the disease is briefly

Return ONLY valid JSON - a single object where keys are the exact trait names and values are the description strings. No markdown, no code fences, no explanation.

Traits:
{trait_list}"""

    try:
        result = subprocess.run(
            ["claude", "-p", "--output-format", "json", "--model", "sonnet", "--max-turns", "1"],
            input=prompt,
            capture_output=True, text=True, timeout=300
        )

        output = result.stdout.strip()

        # Try to extract JSON from the output
        # Claude with --output-format json wraps in {"type":"result","result":"..."}
        try:
            wrapper = json.loads(output)
            if isinstance(wrapper, dict) and "result" in wrapper:
                inner = wrapper["result"]
                # The result might be a JSON string that needs parsing
                if isinstance(inner, str):
                    # Try to find JSON object in the string
                    # Remove any markdown code fences
                    inner = re.sub(r'```json\s*', '', inner)
                    inner = re.sub(r'```\s*', '', inner)
                    descriptions = json.loads(inner)
                elif isinstance(wrapper, dict) and "result" not in wrapper:
                    descriptions = wrapper
                else:
                    descriptions = inner
            elif isinstance(wrapper, dict) and all(isinstance(v, str) for v in wrapper.values()):
                descriptions = wrapper
            else:
                # Try parsing the raw output
                descriptions = json.loads(output)
        except (json.JSONDecodeError, TypeError):
            # Try to find JSON in stdout
            match = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', output, re.DOTALL)
            if match:
                descriptions = json.loads(match.group())
            else:
                print(f"  ERROR: Could not parse response ({len(output)} chars)")
                print(f"  First 200 chars: {output[:200]}")
                continue

        if isinstance(descriptions, dict):
            count = 0
            for trait, desc in descriptions.items():
                if isinstance(desc, str) and len(desc) > 5:
                    existing[trait] = desc
                    count += 1
            print(f"  Got {count} descriptions")
        else:
            print(f"  ERROR: Response is not a dict: {type(descriptions)}")

    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT on batch {batch_idx+1}")
    except Exception as e:
        print(f"  ERROR: {e}")

    # Save after each batch (incremental)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(existing, f, indent=2)
    print(f"  Saved {len(existing)} total descriptions so far")

print(f"\n=== DONE: {len(existing)} descriptions saved to {OUTPUT_FILE} ===")
