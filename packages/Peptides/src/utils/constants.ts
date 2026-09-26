import * as DG from 'datagrok-api/dg';

export enum VIEWER_TYPE {
  SEQUENCE_VARIABILITY_MAP = 'Sequence Variability Map',
  MOST_POTENT_RESIDUES = 'Most Potent Residues',
  LOGO_SUMMARY_TABLE = 'Logo Summary Table',
  DENDROGRAM = 'Dendrogram',
  CLUSTER_MAX_ACTIVITY = 'Active peptide selection',
  MCL = 'MCL',
  SEQUENCE_MUTATION_CLIFFS = 'Sequence Mutation Cliffs',
}

export enum COLUMNS_NAMES {
  SPLIT_COL = '~split',
  ACTIVITY = 'Activity',
  // Deprecated, used for backward compatibility with 1.16.0
  ACTIVITY_SCALED = 'Scaled activity',
  MONOMER = 'AAR',
  POSITION = 'Pos',
  P_VALUE = 'pValue',
  MEAN_DIFFERENCE = 'Mean difference',
  COUNT = 'Count',
  RATIO = 'Ratio',
  MEAN = 'Mean',
  TOTAL_COUNT = '∑'
}

export enum LST_COLUMN_NAMES {
  MEMBERS = 'Members',
  WEB_LOGO = 'WebLogo',
  DISTRIBUTION = 'Distribution',
  MEAN_DIFFERENCE = 'Mean difference',
  P_VALUE = 'P-Value',
  RATIO = 'Ratio',
  CLUSTER = 'Cluster',
}

export enum TAGS {
  MONOMER = 'monomer',
  POSITION = 'pos',
  SEPARATOR = 'separator',
  ALPHABET = 'alphabet',
  FILTER = 'filter',
  SAR_MODE = 'sarMode',
  VISIBLE = 'visible',
  SETTINGS = 'settings',
  CUSTOM_CLUSTER = 'customCluster',
  // UUID = 'pep-uuid',
  MONOMER_POSITION_MODE = 'monomerPositionMode',
  MULTIPLE_VIEWS = 'isMultipleViews',
  IDENTITY_TEMPLATE = 'Identity template',
  SIMILARITY_TEMPLATE = 'Similarity template',
  ANALYSIS_COL = 'isAnalysisCol',
  POSITION_COL = 'isPositionCol',
  M_P_STATS_CACHE = '.MPStatsCache',
  INVARIANT_MAP_COLOR_CACHE = '.InvariantMapColorCache',
  INVARIANT_MAP_COLOR_MAX_CACHE = '.InvariantMapColorMaxCache',
  INVARIANT_MAP_COLOR_MIN_CACHE = '.InvariantMapColorMinCache',
}

export enum SEM_TYPES {
  MONOMER = 'Monomer',
  MACROMOLECULE_DIFFERENCE = 'MacromoleculeDifference',
}

export const COLUMN_NAME = 'ColumnName';

export enum SCALING_METHODS {
  NONE = 'none',
  LG = 'lg',
  MINUS_LG = '-lg',
}

export enum ACTIVITY_TARGET {
  HIGH = 'High',
  LOW = 'Low',
}

export enum SUFFIXES {
  LST = 'lst-', // Logo Summary Table
  MP = 'mp-', // Monomer Position
  MPR = 'mpr-', // Most Potent Residues
  WL = 'wl-', // Web Logo
}

export const AGGREGATION_TYPES = Object.values(DG.AGG)
  .filter((it) => ![DG.AGG.FIRST, DG.AGG.KEY, DG.AGG.PIVOT, DG.AGG.SELECTED_ROWS_COUNT].includes(it));

export const AGG_STATS_MAPPING: {[key: string]: string} = {
  [DG.AGG.TOTAL_COUNT]: 'totalCount',
  [DG.AGG.VALUE_COUNT]: 'valueCount',
  [DG.AGG.UNIQUE_COUNT]: 'uniqueCount',
  [DG.AGG.MISSING_VALUE_COUNT]: 'missingValueCount',
  [DG.AGG.MIN]: 'min',
  [DG.AGG.MAX]: 'max',
  [DG.AGG.SUM]: 'sum',
  [DG.AGG.MED]: 'med',
  [DG.AGG.AVG]: 'avg',
  [DG.AGG.STDEV]: 'stdev',
  [DG.AGG.VARIANCE]: 'variance',
  [DG.AGG.SKEW]: 'skew',
  [DG.AGG.KURT]: 'kurt',
  [DG.AGG.Q1]: 'q1',
  [DG.AGG.Q2]: 'q2',
  [DG.AGG.Q3]: 'q3',
};
