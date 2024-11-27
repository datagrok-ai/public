
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

export enum TrellisAxis {From = 'From', To = 'To'};
export enum TrellisSortByProp {Frequency = 'Frequency', MW = 'MW'};
export enum TrellisSortType {None = 'None', Asc = 'Asc', Desc = 'Desc'};

export const FRAGMENTS_GRID_TOOLTIP = `'Fragment Pairs' grid contains all fragments substitutions for the current molecule in the initial dataset. Molecule pairs with current substitution are shown below in 'Matched Molecular Pairs' grid. Change current row in 'Fragment Pairs' grid to search for molecule pairs with corresponding substitution.`;
export const MATHED_MOLECULAR_PAIRS_TOOLTIP_TRANS = `'Matched Molecular Pairs' grid contains all molecule pairs with current fragment substitution from 'Fragment Pairs' grid. Pinned row is a pair, containing current molecule from initial dataset. Click on any row to select rows with corresponding molecules in initial dataset, selected rows will be pinned on top.`
export const MATHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS = `Click on a non-empty cell in the trellis plot above to see molecule pairs with corresponding fragment substitution. Click on any row to select rows with corresponding molecules in initial dataset, selected rows will be pinned on top.`
export const FRAGMENTS_TAB_TOOLTIP = `Charts in trellis plot cells show mean difference for corresponding activities. Click on any non-empty cell on trellis plot to see molecule pairs with corresponding fragment substitution in 'Matched Molecular Pairs' grid below.`
export const CLIFFS_TAB_TOOLTIP = `2D map where similar molecules are close to each other. Lines connect matched molecular pairs. Arrow points to a molecule with a greater activity value. Use filters on the scatter plot to change activity difference cutoff or switch of/on the lines for any activity. Click on a line to see details in a context panel.`