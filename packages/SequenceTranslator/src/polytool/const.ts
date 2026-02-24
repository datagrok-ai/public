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
export const PT_UI_LINEARIZE_TT = 'Make representation linear if possible';
export const PT_UI_USE_CHIRALITY = 'Chirality engine';
export const PT_UI_HIGHLIGHT_MONOMERS = 'Highlight monomers';
export const PT_UI_DIALOG_CONVERSION = 'Poly Tool Conversion';
export const PT_UI_DIALOG_UNRULE = 'Poly Tool Unrule';
export const PT_UI_DIALOG_ENUMERATION = 'Poly Tool Enumeration';
export const PT_UI_RULES_USED = 'Rules used';

export const PT_ENUM_TYPE_TOOLTIPS: Record<string, string> = {
  'single': 'Each position is enumerated independently. Total results = sum of monomers across all positions.',
  'parallel': 'The i-th result uses the i-th monomer from every position (zip). All positions must have the same number of monomers. Total results = number of monomers per position.',
  'matrix': 'Cartesian product of all positions. Total results = product of monomer counts across all positions.',
  'library': 'Substitutes all monomers from a selected library at a single position.',
};
