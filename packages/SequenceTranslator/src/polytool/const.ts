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
