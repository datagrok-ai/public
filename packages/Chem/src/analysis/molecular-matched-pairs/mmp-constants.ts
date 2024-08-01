import * as DG from 'datagrok-api/dg';

export enum MMP_CONSTRICTIONS {
  CPU = 1E4,
  GPU = 1E5
}

export enum MMP_ERRORS {
  FRAGMENTS_CPU = 'No GPU found and 10,000 molecules is upper limit for MMPa with CPU',
  FRAGMENTS_GPU = 'Upper limit for MMPa with GPU is 100,000 molecules',
  GPU_ABORTED = 'GPU calculations were aborted - faling back to CPU',
  PAIRS = 'Unable to calculate pairs in MMPa analysis',
  GENERATIONS = 'Unable to calculate generations in MMPa analysis'
}

export enum MMP_NAMES {
  FROM = 'From',
  TO = 'To',
  PAIRS = 'Pairs',
  PAIRNUM = '~PairNum',
  PAIRNUM_FROM = '~PairNumFrom',
  PAIRNUM_TO = '~PairNumTo',
  MEANDIFF = 'Mean Difference',
  DIFF = 'Difference',
  VIEW_NAME = 'MMP Analysis',
  TAB_TRANSFORMATIONS = 'Transformations',
  TAB_FRAGMENTS = 'Fragments',
  TAB_CLIFFS = 'Cliffs',
  TAB_GENERATION = 'Generation',
  STRUCT_DIFF_FROM_NAME = '~structDiffFrom',
  STRUCT_DIFF_TO_NAME = '~structDiffTo',
  COLOR = 'color'
}

export const columnsDescriptions: {[key: string]: string} = {
  'Core': 'Common core',
  [MMP_NAMES.FROM]: 'Original fragment',
  [MMP_NAMES.TO]: 'Replacement fragment',
  [MMP_NAMES.PAIRS]: 'Total number of pairs',
  'Generation': 'Transformed molecule',
  'Activity': 'Analyzed Activity/Property',
  'Initial value': 'Initial analyzed activity/property value',
  'Prediction': 'Predicted analyzed activity/property value',
};

export type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

export type MmpInput = {
  table: DG.DataFrame,
  molecules: DG.Column,
  activities: DG.ColumnList,
  fragmentCutoff: number
};
