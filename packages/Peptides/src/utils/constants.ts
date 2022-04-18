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

export enum CATEGORIES {
  OTHER = 'Other',
  ALL = 'All',
}

export enum SEM_TYPES {
  AMINO_ACIDS = 'aminoAcids',
  ALIGNED_SEQUENCE = 'alignedSequence',
}

export const STATS = 'stats';

export const EMBEDDING_STATUS = 'embeddingStatus';

export const PEPTIDES_ANALYSIS = 'isPeptidesAnalysis';

export enum FLAGS {
  CELL_CHANGING = 'isCellChanging',
}
