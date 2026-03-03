#name: Number Antibody Sequences (ANARCI)
#description: [Legacy] Assigns antibody numbering (IMGT/Kabat/Chothia/AHo) using ANARCI directly
#language: python
#environment: channels: [conda-forge, bioconda, defaults], dependencies: [python=3.9, pip, hmmer=3.3.2, anarci]
#input: dataframe df
#input: column seqCol {semType: Macromolecule}
#input: string scheme {choices: ["imgt", "kabat", "chothia", "aho"]} [Numbering scheme]
#output: dataframe result

import json
import pandas as pd

try:
    from anarci import anarci, number
except ImportError:
    raise ImportError(
        'ANARCI is not installed. Please install it with:\n'
        '  conda install -c bioconda anarci hmmer\n'
        'or:\n'
        '  pip install anarci\n'
        'HMMER must also be available on PATH.'
    )

# IMGT region boundaries (used for all schemes after numbering)
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
    # If it looks like FASTA single-letter (only AA chars and gaps)
    s = s.replace('-', '').replace('.', '')
    # Filter to valid amino acid characters
    valid = set('ACDEFGHIKLMNPQRSTVWY')
    return ''.join(c for c in s.upper() if c in valid)


def number_sequences(sequences, scheme_name):
    """Number a list of sequences using ANARCI."""
    # Prepare input for ANARCI
    seq_list = [(f'seq_{i}', seq) for i, seq in enumerate(sequences) if seq]

    if not seq_list:
        return [], [], []

    # Run ANARCI
    numbered_seqs, alignment_details, hit_tables = anarci(
        seq_list,
        scheme=scheme_name,
        output=False,
        allow=set(['H', 'K', 'L']),
    )

    return numbered_seqs, alignment_details, hit_tables


# Extract sequences from the input column
col_name = seqCol
sequences = []
for i in range(df.shape[0]):
    raw = df[col_name].iloc[i]
    sequences.append(extract_sequence(str(raw) if raw is not None else ''))

# Number all sequences
numbered_seqs, alignment_details, hit_tables = number_sequences(sequences, scheme)

# Process results
position_names_list = []
chain_types = []
region_annotations_list = []
numbering_results = []

regions_def = SCHEME_REGIONS.get(scheme, IMGT_REGIONS)

for idx in range(len(sequences)):
    if numbered_seqs[idx] is None or len(numbered_seqs[idx]) == 0:
        position_names_list.append('')
        chain_types.append('')
        region_annotations_list.append('[]')
        numbering_results.append('')
        continue

    # Get the best hit (first domain)
    numbering = numbered_seqs[idx][0]
    chain_info = alignment_details[idx][0]

    if numbering is None or chain_info is None:
        position_names_list.append('')
        chain_types.append('')
        region_annotations_list.append('[]')
        numbering_results.append('')
        continue

    # Extract chain type
    chain_type_code = chain_info['chain_type']
    if chain_type_code == 'H':
        chain_type = 'Heavy'
    elif chain_type_code in ('K', 'L'):
        chain_type = 'Light'
    else:
        chain_type = 'Unknown'

    # Build position names from numbering
    pos_names = []
    for (pos_num, insertion_code), aa in numbering:
        if aa != '-':
            pos_name = str(pos_num)
            if insertion_code != ' ':
                pos_name += insertion_code
            pos_names.append(pos_name)

    position_names_list.append(', '.join(pos_names))
    chain_types.append(chain_type)

    # Build region annotations
    chain_key = chain_type if chain_type in regions_def else 'Heavy'
    region_defs = regions_def.get(chain_key, [])
    annotations = []
    for region_name, start, end in region_defs:
        region_type = 'CDR' if 'CDR' in region_name else 'FR'
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

    # Full numbering as JSON for debugging/reference
    numbering_json = json.dumps([
        {'position': f'{pos_num}{insertion_code.strip()}', 'aa': aa}
        for (pos_num, insertion_code), aa in numbering
        if aa != '-'
    ])
    numbering_results.append(numbering_json)

# Build result DataFrame
result = pd.DataFrame({
    'position_names': position_names_list,
    'chain_type': chain_types,
    'annotations_json': region_annotations_list,
    'numbering_detail': numbering_results,
})
