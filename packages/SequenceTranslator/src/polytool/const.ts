import {HELM_REQUIRED_FIELD} from '@datagrok-libraries/bio/src/utils/const';

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
  [HELM_REQUIRED_FIELD.SYMBOL]: 'Short Name',
  [HELM_REQUIRED_FIELD.NAME]: 'Medium Name',
  [HELM_REQUIRED_FIELD.SMILES]: 'SMILES',
};

export const R_GROUP_BLOCK_DUMMY = [
  {
    'capGroupSMILES': '[*:1][H]',
    'alternateId': 'R1-H',
    'capGroupName': 'H',
    'label': 'R1'
  },
  {
    'capGroupSMILES': 'O[*:2]',
    'alternateId': 'R2-OH',
    'capGroupName': 'OH',
    'label': 'R2'
  },
  {
    'capGroupSMILES': '[*:3][H]',
    'alternateId': 'R3-H',
    'capGroupName': 'H',
    'label': 'R3'
  }
];
