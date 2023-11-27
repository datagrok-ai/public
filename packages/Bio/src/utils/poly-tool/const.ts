import {HELM_REQUIRED_FIELDS} from '@datagrok-libraries/bio/src/utils/const';

export const enum HELM_WRAPPER {
  LEFT = 'PEPTIDE1{',
  RIGHT = '}$$$$',
}
export const ALL_MONOMERS = '<All>';

export const enum TRANSFORMATION_TYPE {
  CYCLIZATION = 'Cyclization',
}

export const enum CYCLIZATION_TYPE {
  NO = 'N-O',
  NCys = 'N-Cys',
  R3 = 'R3-R3',
}

export const helmFieldsToPolyToolInputFields = {
  [HELM_REQUIRED_FIELDS.SYMBOL]: 'Short Name',
  [HELM_REQUIRED_FIELDS.NAME]: 'Medium Name',
  [HELM_REQUIRED_FIELDS.SMILES]: 'SMILES',
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
