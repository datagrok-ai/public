export const MODIFICATIONS: {[index: string]: {molecularWeight: number, left: string, right: string}} = {
  '(invabasic)': {
    molecularWeight: 118.13,
    left: 'O[C@@H]1C[C@@H]O[C@H]1CO',
    right: 'O[C@@H]1C[C@@H]O[C@H]1CO',
  },
  '(GalNAc-2-JNJ)': {
    molecularWeight: 1273.3,
    left: 'C(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)NC(=O)CCCC(=O)NCC(O)CO',
    right: 'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)'+
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)',
  },
};

export const standardPhosphateLinkSmiles = 'OP(=O)(O)O';
