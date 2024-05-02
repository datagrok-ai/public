export const MMP_COLNAME_FROM = 'From';
export const MMP_COLNAME_TO = 'To';
export const MMP_COLNAME_PAIRS = 'Pairs';
export const MMP_COL_PAIRNUM = '~PairNum';
export const MMP_COL_PAIRNUM_FROM = '~PairNumFrom';
export const MMP_COL_PAIRNUM_TO = '~PairNumTo';
export const MMP_COLNAME_MEANDIFF = 'Mean Difference';
export const MMP_COLNAME_DIFF = 'Difference';
export const MMP_VIEW_NAME = 'MMP Analysis';
export const MMP_TAB_TRANSFORMATIONS = 'Transformations';
export const MMP_TAB_FRAGMENTS = 'Fragments';
export const MMP_TAB_CLIFFS = 'Cliffs';
export const MMP_TAB_GENERATION = 'Generation';
export const MMP_STRUCT_DIFF_FROM_NAME = '~structDiffFrom';
export const MMP_STRUCT_DIFF_TO_NAME = '~structDiffTo';
export const MMP_COLOR = 'color';

export const columnsDescriptions: {[key: string]: string}  = {
    "Core": "Common core",
    [MMP_COLNAME_FROM]: "Original fragment",
    [MMP_COLNAME_TO]: "Replacement fragment",
    [MMP_COLNAME_PAIRS]: 'Total number of pairs',
    "Generation": "Transformed molecule",
    "Activity": "Analyzed Activity/Property",
    "Initial value": "Initial analyzed activity/property value",
    "Prediction": "Predicted analyzed activity/property value",
  }
