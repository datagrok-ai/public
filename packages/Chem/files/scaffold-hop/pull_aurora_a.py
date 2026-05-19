"""Pull all ChEMBL Aurora kinase A (CHEMBL4722) IC50 activities, deduplicate
by compound, mask ~25% of pIC50 to NA, and write a CSV ready for the MMP
imputation + scaffold-hopping workflow.

Resumable: each page is cached in ``.chembl_cache/page_<offset>.json``. Re-running
the script skips already-cached offsets and only fetches gaps. Useful because the
ChEMBL API is flaky from this network (intermittent HTTP 500s, timeouts).

Uses stdlib only — no third-party deps.

Output columns:
  smiles      — canonical SMILES from ChEMBL
  CHEMBL_ID   — compound ChEMBL ID (for traceability)
  pIC50       — median pIC50 across all measurements for this compound;
                "NA" for the masked subset
  pIC50_true  — same as pIC50 BUT with the masked rows showing the true
                values; used by the imputation holdout test as the
                ground-truth column

Usage:
  python pull_aurora_a.py                  # default: fetch + write CSV
  python pull_aurora_a.py --write-only     # only write CSV from existing cache (no fetch)
"""
import csv
import json
import math
import os
import random
import statistics
import subprocess
import sys
import time

TARGET_ID = 'CHEMBL4722'                           # Aurora kinase A
STANDARD_TYPE = 'IC50'
PAGE_SIZE = 200                                    # 1000 timed out; 200 is reliable
SLEEP_BETWEEN_PAGES = 0.3                          # be nice to ChEMBL
MAX_RETRIES = 4                                    # per page
RANDOM_SEED = 42                                   # deterministic mask
MASK_FRACTION = 0.25                               # ~25% of compounds get NA

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(OUT_DIR, '.chembl_cache')
OUT_CSV = os.path.join(OUT_DIR, 'aurora_a_chembl.csv')


def fetch_page(offset, page_size=PAGE_SIZE):
    """Fetch one ChEMBL page via curl. Returns parsed JSON or None on failure."""
    url = (f'https://www.ebi.ac.uk/chembl/api/data/activity.json'
           f'?target_chembl_id={TARGET_ID}'
           f'&standard_type={STANDARD_TYPE}'
           f'&limit={page_size}'
           f'&offset={offset}')
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            raw = subprocess.check_output(
                ['curl', '-sS', '--fail', '--max-time', '120', url],
                timeout=150, stderr=subprocess.STDOUT,
            )
            return json.loads(raw.decode('utf-8'))
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
                json.JSONDecodeError) as e:
            last_err = e
            sys.stderr.write(f'    attempt {attempt}/{MAX_RETRIES} failed: '
                             f'{type(e).__name__}; retrying...\n')
            time.sleep(1.0 + attempt)
    sys.stderr.write(f'    GAVE UP on offset={offset}: {last_err}\n')
    return None


def load_cached_pages():
    """Return list of (offset, parsed_json) for every cache file on disk."""
    pages = []
    if not os.path.isdir(CACHE_DIR):
        return pages
    for fname in sorted(os.listdir(CACHE_DIR)):
        if not fname.endswith('.json'):
            continue
        path = os.path.join(CACHE_DIR, fname)
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            sys.stderr.write(f'  skipping bad cache file {fname}: {e}\n')
            continue
        pm = data.get('page_meta', {})
        off = pm.get('offset')
        if off is None:
            continue
        pages.append((off, fname, data))
    return pages


def fetch_all_with_cache():
    """Fill in gaps in the cache, saving each new page to disk immediately."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cached = load_cached_pages()
    sys.stderr.write(f'  cached pages: {len(cached)}\n')

    # Discover total_count from any existing page; fall back to a probe.
    total = None
    for off, _fname, data in cached:
        total = data.get('page_meta', {}).get('total_count')
        if total:
            break
    if total is None:
        sys.stderr.write('  no cache → probing for total_count\n')
        probe = fetch_page(0, page_size=1)
        if probe is None:
            sys.stderr.write('  PROBE FAILED — cannot discover total_count\n')
            return cached
        total = probe.get('page_meta', {}).get('total_count', 0)
        sys.stderr.write(f'  total_count={total}\n')
    sys.stderr.write(f'  ChEMBL reports total_count={total}\n')

    # Build set of offsets already covered (we treat each cached file as
    # covering [offset, offset+limit). Standard pagination at limit=PAGE_SIZE
    # then asks for any offset whose window is not yet covered.)
    covered = set()
    for off, _fname, data in cached:
        lim = data.get('page_meta', {}).get('limit', PAGE_SIZE)
        for o in range(off, off + lim):
            covered.add(o)

    # Standard pagination: offsets 0, PAGE_SIZE, 2*PAGE_SIZE, …
    # A page-offset is "needed" if any of its records [off, off+PAGE_SIZE) is
    # uncovered.
    missing_offsets = []
    for off in range(0, total, PAGE_SIZE):
        if any(o not in covered for o in range(off, off + PAGE_SIZE)):
            missing_offsets.append(off)

    if not missing_offsets:
        sys.stderr.write('  cache already complete\n')
        return cached
    sys.stderr.write(f'  missing {len(missing_offsets)} pages: '
                     f'{missing_offsets[:5]}…{missing_offsets[-3:]}\n')

    for i, off in enumerate(missing_offsets, 1):
        sys.stderr.write(f'  [{i}/{len(missing_offsets)}] fetching offset={off}... ')
        sys.stderr.flush()
        page = fetch_page(off, page_size=PAGE_SIZE)
        if page is None:
            sys.stderr.write(f'    skip\n')
            continue
        acts = page.get('activities', [])
        sys.stderr.write(f'got {len(acts)} records\n')
        path = os.path.join(CACHE_DIR, f'page_{off}.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(page, f)
        cached.append((off, f'page_{off}.json', page))
        time.sleep(SLEEP_BETWEEN_PAGES)

    return cached


def to_pic50(record):
    """Convert one activity record to pIC50. Prefers ChEMBL's pre-computed
    pchembl_value when available; falls back to IC50 conversion. Returns
    None on invalid records (missing value, non-numeric units, etc.)."""
    # ChEMBL pre-computes pchembl_value for high-confidence numeric records.
    pchembl = record.get('pchembl_value')
    if pchembl:
        try:
            return float(pchembl)
        except (TypeError, ValueError):
            pass
    # Fall back: compute from standard_value + standard_units.
    val = record.get('standard_value')
    units = record.get('standard_units')
    if val is None or units is None:
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    if v <= 0:
        return None
    # Normalise to molar concentration.
    unit_to_molar = {
        'nM': 1e-9, 'nm': 1e-9,
        'uM': 1e-6, 'um': 1e-6, 'µM': 1e-6,
        'mM': 1e-3,
        'pM': 1e-12,
        'M': 1.0,
    }
    factor = unit_to_molar.get(units)
    if factor is None:
        return None
    molar = v * factor
    return -math.log10(molar)


def aggregate_records(pages):
    """pages = list of (offset, fname, parsed_json) → flat list of activity dicts."""
    records = []
    for _off, _fname, data in pages:
        records.extend(data.get('activities', []))
    return records


def main():
    write_only = '--write-only' in sys.argv
    sys.stderr.write(f'Pulling {STANDARD_TYPE} activities for {TARGET_ID}\n')

    if write_only:
        pages = load_cached_pages()
        sys.stderr.write(f'  --write-only: using {len(pages)} cached pages\n')
    else:
        pages = fetch_all_with_cache()

    records = aggregate_records(pages)
    sys.stderr.write(f'\nTotal raw records (incl. duplicates across pages): {len(records)}\n')
    if not records:
        sys.stderr.write('No records — aborting.\n')
        sys.exit(1)

    # Group by compound. Keep median pIC50 per (compound, smiles) pair.
    per_compound = {}  # chembl_id -> {'smiles': ..., 'pic50s': [...]}
    skipped_no_smiles = 0
    skipped_no_pic50 = 0
    for rec in records:
        cid = rec.get('molecule_chembl_id')
        smi = rec.get('canonical_smiles')
        if not cid:
            continue
        if not smi:
            skipped_no_smiles += 1
            continue
        pic50 = to_pic50(rec)
        if pic50 is None:
            skipped_no_pic50 += 1
            continue
        if cid not in per_compound:
            per_compound[cid] = {'smiles': smi, 'pic50s': []}
        per_compound[cid]['pic50s'].append(pic50)

    # Collapse: take median pIC50 per compound.
    rows = []
    for cid, info in per_compound.items():
        if not info['pic50s']:
            continue
        rows.append({
            'smiles': info['smiles'],
            'CHEMBL_ID': cid,
            'pIC50': statistics.median(info['pic50s']),
        })

    sys.stderr.write(f'\nSkipped:\n')
    sys.stderr.write(f'  no SMILES:        {skipped_no_smiles}\n')
    sys.stderr.write(f'  no/invalid pIC50: {skipped_no_pic50}\n')
    sys.stderr.write(f'\nUnique compounds after dedup: {len(rows)}\n')
    if rows:
        pic50s = [r['pIC50'] for r in rows]
        sys.stderr.write(
            f'pIC50 range: {min(pic50s):.2f}–{max(pic50s):.2f}, '
            f'median={statistics.median(pic50s):.2f}, '
            f'stdev={statistics.stdev(pic50s):.2f}\n')

    # Deterministic ~25% mask. The masked rows show NA in the pIC50 column
    # AND the true value in pIC50_true — the imputation holdout test reads
    # pIC50_true to compute MAE/R² on the masked subset.
    random.seed(RANDOM_SEED)
    rows.sort(key=lambda r: r['CHEMBL_ID'])           # deterministic ordering
    n_masked = int(len(rows) * MASK_FRACTION)
    mask_idxs = set(random.sample(range(len(rows)), n_masked))
    for i, r in enumerate(rows):
        r['pIC50_true'] = r['pIC50']
        if i in mask_idxs:
            r['pIC50'] = 'NA'
    sys.stderr.write(f'\nMasked {n_masked} of {len(rows)} compounds ({MASK_FRACTION*100:.0f}%)\n')

    # Write CSV.
    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['smiles', 'CHEMBL_ID', 'pIC50', 'pIC50_true'])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    sys.stderr.write(f'\nWrote {OUT_CSV}\n')


if __name__ == '__main__':
    main()
