
export enum MMP_NAMES {
  FROM = 'From',
  TO = 'To',
  PAIRS_COUNT = 'Total count',
  PAIRNUM = '~PairNum',
  PAIR_SORT = '~PairSort',
  PAIRNUM_FROM = '~PairNumFrom',
  PAIRNUM_TO = '~PairNumTo',
  MEANDIFF = '\u0394',
  DIFF = '\u0394',
  VIEW_NAME = 'MMP Analysis',
  TAB_TRANSFORMATIONS = 'Substitutions',
  TAB_FRAGMENTS = 'Fragments',
  TAB_CLIFFS = 'Cliffs',
  TAB_GENERATION = 'Generation',
  COLOR = 'color',
  SMI1 = '~smi1',
  SMI2 = '~smi2',
  RULENUM = '~ruleNum',
  FRAGMENTS_GRID = 'Fragments',
  PAIRS_GRID = 'Molecule Pairs',
  STRUCTURE = 'Structure',
  INITIAL_VALUE = 'Initial value',
  ACTIVITY = 'Activity',
  CORE = 'Core',
  CORE_NUM = '~Core',
  PREDICTION = 'Prediction',
  GENERATIONS = 'Generation',
  NEW_MOLECULE = 'Predicted molecule',
  OBSERVED = 'Observed',
  PREDICTED = 'Predicted',
  PROPERTY_TYPE = 'Property type',
  DELTA_ACTIVITY = '\u0394 activity',
}

export enum SHOW_FRAGS_MODE {
  All = 'All',
  Current = 'Current molecule'
}

export const FRAGMENTS_GRID_HEADER_DESCRIPTIONS: {[key: string]: string} = {
  [MMP_NAMES.FROM]: 'Initial fragment',
  [MMP_NAMES.TO]: 'Fragment after substitution',
  [MMP_NAMES.PAIRS_COUNT]: 'Total number of pairs for substitution identified in a dataset',
};

export const PAIRS_GRID_HEADER_DESCRIPTIONS: {[key: string]: string} = {
  [MMP_NAMES.FROM]: 'Initial molecule',
  [MMP_NAMES.TO]: 'Molecule after substitution',
};

export const GENERATIONS_GRID_HEADER_DESCRIPTIONS: {[key: string]: string} = {
  [MMP_NAMES.STRUCTURE]: 'Structure from initial dataset',
  [MMP_NAMES.CORE]: 'Conservative fragment',
  [MMP_NAMES.FROM]: 'Initial fragment',
  [MMP_NAMES.TO]: 'Fragment after substitution',
  [MMP_NAMES.NEW_MOLECULE]: 'New molecule obtained with fragments substitution',
  [MMP_NAMES.PROPERTY_TYPE]: 'Property type',
  [MMP_NAMES.OBSERVED]: 'Observed activity value',
  [MMP_NAMES.PREDICTED]: 'Predicted activity value based on MMP rules (initial activity + mean difference)',
  [MMP_NAMES.DELTA_ACTIVITY]: 'Activity difference',
  [MMP_NAMES.PREDICTION]: 'True if molecule does not exists in the initial dataset',
};

export enum TrellisAxis {From = 'From', To = 'To'};
export enum TrellisSortByProp {Frequency = 'Frequency', MW = 'MW', None = 'None'};
export enum TrellisSortType {Asc = 'Asc', Desc = 'Desc'};

export const FRAGMENTS_GRID_TOOLTIP = `'Fragment Pairs' grid contains all fragment substitutions found in the dataset.`;
export const MATCHED_MOLECULAR_PAIRS_TOOLTIP_TRANS = `'Molecular Pairs' grid shows pairs of molecules with found substitutions. Click any row to pin corresponding molecules on top of initial dataset.`;
export const MATCHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS = `Click a non-empty cell in the trellis plot above to filter molecule pairs with corresponding fragment substitution.`;
export const MATCHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS = `Click any row to zoom to corresponding pair on scatter plot`;
export const FRAGMENTS_TAB_TOOLTIP = `Charts in trellis plot cells show mean difference for corresponding activities. Click any non-empty cell on trellis plot to filter molecule pairs with corresponding fragment substitution in 'Matched Molecular Pairs' grid below.`;
export const CLIFFS_TAB_TOOLTIP = `2D map where similar molecules are close to each other. Lines connect matched molecular pairs. Arrow points to a molecule with a greater activity value. Use filters on the scatter plot to change activity difference cutoff or switch of/on the lines for any activity. Click on a line to see molecule pair details in a context panel.`;
export const GENERATIONS_TAB_TOOLTIP = `Generation of new molecules based on obtained rules. For each molecule from initial dataset and for each selected activity the most optimal substitution is shown with associated prediction results for all activities. We specify if initial or generated molecule is optimal. Open context panel to estimate generation model predictive ability on observed vs predicted plot.`;

export const MMP_CONTEXT_PANE_CLASS = 'mmp-context-pane';

export const MOL_CANVAS_WIDTH = 170;
export const MOL_CANVAS_HEIGHT = 75;

export enum SCALING_METHODS {
  NONE = 'none',
  LG = 'lg',
  MINUS_LG = '-lg',
}
