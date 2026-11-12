export enum COLUMNS_NAMES {
  SPLIT_COL = '~split',
  ACTIVITY = '~activity',
  ACTIVITY_SCALED = 'activity_scaled',
  ALIGNED_SEQUENCE = '~aligned_sequence',
  AMINO_ACID_RESIDUE = 'AAR',
  POSITION = 'Pos',
  P_VALUE = 'pValue',
  MEAN_DIFFERENCE = 'Mean difference',
}

export enum TAGS {
  AAR = 'AAR',
  POSITION = 'Pos',
  SELECTION = 'selection',
}

export enum SEM_TYPES {
  MONOMER = 'Monomer',
  MACROMOLECULE_DIFFERENCE = 'MacromoleculeDifference',
  ACTIVITY = 'activity',
  ACTIVITY_SCALED = 'activityScaled'
}

export const MSA = 'MSA';

export const STATS = 'stats';

export enum FLAGS {
  CELL_CHANGING = 'isCellChanging',
}

export const groupDescription: { [key: string]: { 'description': string, aminoAcids: string[] } } = {
  'PC': {'description': 'Positive Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['R', 'H', 'K']},
  'NC': {'description': 'Negative Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['D', 'E']},
  'U': {'description': 'Amino Acids with Polar Uncharged Side Chains', 'aminoAcids': ['S', 'T', 'N', 'Q']},
  'SC': {'description': 'Special Cases', 'aminoAcids': ['C', 'U', 'G', 'P']},
  'H': {
    'description': 'Amino Acids with Hydrophobic Side Chain',
    'aminoAcids': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
  },
  '-': {'description': 'Unknown Amino Acid', 'aminoAcids': ['-']},
};

export namespace PEPSEA {
  export const SEPARATOR = '.';
}

export const kalignVersion = '3.3.1';

/** meta.role value for dynamically discovered MSA engine functions */
export const MSA_ENGINE_ROLE = 'sequenceMSA';
