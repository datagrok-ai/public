import {HELM_REQUIRED_FIELD} from '@datagrok-libraries/bio/src/utils/const';

export const helmFieldsToPolyToolInputFields = {
  [HELM_REQUIRED_FIELD.SYMBOL]: 'Short Name',
  [HELM_REQUIRED_FIELD.NAME]: 'Medium Name',
  [HELM_REQUIRED_FIELD.SMILES]: 'SMILES',
};

export const R_GROUP_BLOCK_DUMMY = [
  {
    'capGroupSmiles': '[*:1][H]',
    'alternateId': 'R1-H',
    'capGroupName': 'H',
    'label': 'R1'
  },
  {
    'capGroupSmiles': 'O[*:2]',
    'alternateId': 'R2-OH',
    'capGroupName': 'OH',
    'label': 'R2'
  },
  {
    'capGroupSmiles': '[*:3][H]',
    'alternateId': 'R3-H',
    'capGroupName': 'H',
    'label': 'R3'
  }
];

export const PT_ERROR_DATAFRAME = 'No dataframe with macromolecule columns open';
export const PT_WARNING_COLUMN = 'No marcomolecule column chosen!';

export const PT_UI_GET_HELM = 'Get HELM';
export const PT_UI_ADD_HELM = 'Add HELM column';
export const PT_UI_LINEARIZE = 'Linearize';
export const PT_UI_LINEARIZE_TT = 'Draw sequences without cycles or branches as a straight chain, ' +
  'the others as a 2D depiction';
export const PT_UI_HIGHLIGHT_MONOMERS = 'Highlight monomers';
export const PT_UI_HIGHLIGHT_MONOMERS_TT = 'Color each molecule\'s atoms by the monomer they come from ' +
  '(the monomer\'s color in the library)';
export const PT_UI_COLUMN_TT = 'The sequence column (custom notation) to convert';
export const PT_UI_KEEP_ORIGINAL_TT = 'Add the original sequence as the first result';
export const PT_UI_TO_ATOMIC_LEVEL_TT = 'Add a molecule column for the enumerated sequences';
export const PT_UI_TRIVIAL_NAME_TT = 'An ID column for the results: this name and each result\'s changes, ' +
  'e.g. Pep1-A3G (old monomer, position, new monomer)';
export const PT_UI_APPEND_TO_TABLE_TT = 'Add the results to this table instead of a new one';
export const PT_UI_DIALOG_CONVERSION = 'PolyTool Conversion';
export const PT_UI_DIALOG_UNRULE = 'PolyTool Unrule';
export const PT_HELM_UI_DIALOG_ENUMERATION = 'PolyTool Helm Enumeration';
export const PT_UI_RULES_USED = 'Rules used';

export const PT_ENUM_TYPE_TOOLTIPS: Record<string, string> = {
  'single': 'Each position is enumerated independently. Total results = sum of monomers across all positions.',
  'parallel': 'The i-th result uses the i-th monomer from every position (zip). All positions must have the same number of monomers. Total results = number of monomers per position.',
  'matrix': 'Cartesian product of all positions. Total results = product of monomer counts across all positions.',
};
