
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
  COLOR = 'color',
  SMI1 = '~smi1',
  SMI2 = '~smi2',
  RULENUM = '~ruleNum'
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
