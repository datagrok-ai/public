export enum COLUMNS_NAMES {
  SPLIT_COL = '~split',
  ACTIVITY_SCALED = 'activity_scaled',
  MONOMER = 'AAR',
  POSITION = 'Pos',
  P_VALUE = 'pValue',
  MEAN_DIFFERENCE = 'Mean difference',
  COUNT = 'Count',
  RATIO = 'Ratio',
  CLUSTERS = 'clusters',
  MEMBERS = 'Members',
}

export enum CATEGORIES {
  OTHER = 'Other',
  ALL = 'All',
}

export enum TAGS {
  MONOMER = 'monomer',
  POSITION = 'pos',
  SEPARATOR = 'separator',
  SELECTION = 'selection',
  ALPHABET = 'alphabet',
  FILTER = 'filter',
  CLUSTERS = 'clusters',
  SAR_MODE = 'sarMode',
  CLUSTER_SELECTION = 'clusterSelection',
  VISIBLE = 'visible',
  SETTINGS = 'settings',
}

export enum SEM_TYPES {
  MONOMER = 'Monomer',
  MACROMOLECULE_DIFFERENCE = 'MacromoleculeDifference',
}

export const STATS = 'stats';

export const EMBEDDING_STATUS = 'embeddingStatus';

export const PEPTIDES_ANALYSIS = 'isPeptidesAnalysis';

export enum FLAGS {
  CELL_CHANGING = 'isCellChanging',
}

export const aarGroups = {
  'R': 'PC', 'H': 'PC', 'K': 'PC',
  'D': 'NC', 'E': 'NC',
  'S': 'U', 'T': 'U', 'N': 'U', 'Q': 'U',
  'C': 'SC', 'U': 'SC', 'G': 'SC', 'P': 'SC',
  'A': 'H', 'V': 'H', 'I': 'H', 'L': 'H', 'M': 'H', 'F': 'H', 'Y': 'H', 'W': 'H',
  '-': '-',
};

export const groupDescription: {[key: string]: {'description': string, aminoAcids: string[]}} = {
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
