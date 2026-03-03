#name: Number Antibody Sequences (AbNumber)
#description: [Legacy] Assigns antibody numbering (IMGT/Kabat/Chothia/AHo) using AbNumber
#language: python
#environment: channels: [conda-forge, bioconda, defaults], dependencies: [python=3.9, pip, hmmer=3.3.2, pip: [abnumber]]
#input: dataframe df
#input: column seqCol {semType: Macromolecule}
#input: string scheme {choices: ["imgt", "kabat", "chothia", "aho"]} [Numbering scheme]
#output: dataframe result

import json
import pandas as pd

try:
    from abnumber import Chain
except ImportError:
    raise ImportError(
        'AbNumber is not installed. Please install it with:\n'
        '  conda install -c bioconda abnumber\n'
        'or:\n'
        '  pip install abnumber\n'
        'HMMER must also be available on PATH.'
    )

# Region definitions per scheme, keyed by (scheme, chain_group).
# chain_group is 'Heavy' or 'Light'.
# Each region: (name, start_position, end_position).
IMGT_REGIONS = {
    'Heavy': [
        ('FR1', 1, 26), ('CDR1', 27, 38), ('FR2', 39, 55),
        ('CDR2', 56, 65), ('FR3', 66, 104), ('CDR3', 105, 117), ('FR4', 118, 128),
    ],
    'Light': [
        ('FR1', 1, 26), ('CDR1', 27, 38), ('FR2', 39, 55),
        ('CDR2', 56, 65), ('FR3', 66, 104), ('CDR3', 105, 117), ('FR4', 118, 127),
    ],
}

KABAT_REGIONS = {
    'Heavy': [
        ('FR1', 1, 30), ('CDR1', 31, 35), ('FR2', 36, 49),
        ('CDR2', 50, 65), ('FR3', 66, 94), ('CDR3', 95, 102), ('FR4', 103, 113),
    ],
    'Light': [
        ('FR1', 1, 23), ('CDR1', 24, 34), ('FR2', 35, 49),
        ('CDR2', 50, 56), ('FR3', 57, 88), ('CDR3', 89, 97), ('FR4', 98, 107),
    ],
}

CHOTHIA_REGIONS = {
    'Heavy': [
        ('FR1', 1, 25), ('CDR1', 26, 32), ('FR2', 33, 51),
        ('CDR2', 52, 56), ('FR3', 57, 94), ('CDR3', 95, 102), ('FR4', 103, 113),
    ],
    'Light': [
        ('FR1', 1, 25), ('CDR1', 26, 32), ('FR2', 33, 49),
        ('CDR2', 50, 52), ('FR3', 53, 90), ('CDR3', 91, 96), ('FR4', 97, 107),
    ],
}

AHO_REGIONS = {
    'Heavy': [
        ('FR1', 1, 24), ('CDR1', 25, 40), ('FR2', 41, 55),
        ('CDR2', 56, 78), ('FR3', 79, 108), ('CDR3', 109, 138), ('FR4', 139, 149),
    ],
    'Light': [
        ('FR1', 1, 24), ('CDR1', 25, 40), ('FR2', 41, 55),
        ('CDR2', 56, 78), ('FR3', 79, 108), ('CDR3', 109, 138), ('FR4', 139, 149),
    ],
}

SCHEME_REGIONS = {
    'imgt': IMGT_REGIONS,
    'kabat': KABAT_REGIONS,
    'chothia': CHOTHIA_REGIONS,
    'aho': AHO_REGIONS,
}


def extract_sequence(raw_seq):
    """Extract single-letter amino acid sequence from various formats."""
    if not raw_seq or not isinstance(raw_seq, str):
        return ''
    s = raw_seq.strip()
    s = s.replace('-', '').replace('.', '')
    valid = set('ACDEFGHIKLMNPQRSTVWY')
    return ''.join(c for c in s.upper() if c in valid)


def number_single_sequence(seq, scheme_name):
    """Number a single sequence using AbNumber's Chain object.

    Returns (chain, error_msg) where chain is an abnumber.Chain or None.
    """
    if not seq:
        return None, 'empty sequence'
    try:
        chain = Chain(seq, scheme=scheme_name)
        return chain, None
    except Exception as e:
        return None, str(e)


# Extract sequences from the input column
col_name = seqCol
sequences = []
for i in range(df.shape[0]):
    raw = df[col_name].iloc[i]
    sequences.append(extract_sequence(str(raw) if raw is not None else ''))

# Process each sequence
position_names_list = []
chain_types = []
region_annotations_list = []
numbering_results = []

regions_def = SCHEME_REGIONS.get(scheme, IMGT_REGIONS)

for seq in sequences:
    chain, err = number_single_sequence(seq, scheme)

    if chain is None:
        position_names_list.append('')
        chain_types.append('')
        region_annotations_list.append('[]')
        numbering_results.append('')
        continue

    # Determine chain type
    ct = chain.chain_type  # 'H', 'K', or 'L'
    if ct == 'H':
        chain_type = 'Heavy'
    elif ct in ('K', 'L'):
        chain_type = 'Light'
    else:
        chain_type = 'Unknown'

    # Build position names from numbering
    pos_names = []
    numbering_detail = []
    for position, aa in chain:
        pos_name = str(position.number)
        if position.letter and position.letter.strip():
            pos_name += position.letter.strip()
        pos_names.append(pos_name)
        numbering_detail.append({'position': pos_name, 'aa': aa})

    position_names_list.append(', '.join(pos_names))

    chain_types.append(chain_type)

    # Build region annotations
    chain_key = chain_type if chain_type in regions_def else 'Heavy'
    region_defs = regions_def.get(chain_key, [])
    annotations = []
    for region_name, start, end in region_defs:
        annotations.append({
            'id': f'{scheme}-{chain_type}-{region_name}'.lower(),
            'name': region_name,
            'description': f'{region_name} ({scheme.upper()} {start}-{end})',
            'start': str(start),
            'end': str(end),
            'visualType': 'region',
            'category': 'structure',
            'sourceScheme': scheme.upper(),
            'autoGenerated': True,
        })

    region_annotations_list.append(json.dumps(annotations))
    numbering_results.append(json.dumps(numbering_detail))

# Build result DataFrame
result = pd.DataFrame({
    'position_names': position_names_list,
    'chain_type': chain_types,
    'annotations_json': region_annotations_list,
    'numbering_detail': numbering_results,
})
