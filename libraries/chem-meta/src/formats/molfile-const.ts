export const R_GROUP_ELEMENT_SYMBOL = 'R#';
export const HYDROGEN_SYMBOL = 'H';

export enum V2K_CONST {
  TYPE = 'V2000',
  NUM_OF_HEADER_LINES = 3,

  NUM_OF_COUNTS_DIGITS = 3,
  ATOM_TYPE_COL = 4,
  FIRST_BONDED_ATOM_COL = 1,
  BOND_TYPE_COL = 3,

  RGP_SHIFT = 8,
  MAX_ATOM_COUNT = 999,

  RGP_LINE_START = 'M  RGP',
  ATOM_ALIAS_LINE_START = 'A  ',

  END = 'M  END',
}

export const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;

export const enum V3K_CONST {
  TYPE = 'V3000',
  HEADER_SECOND_LINE = '  0  0  0  0  0  0            999 V3000\n',

  BEGIN_DATA_LINE = 'M  V30 ',

  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  /** Index shift from the beginning of the 'COUNTS' line to the number of atoms  */
  COUNTS_SHIFT = 14,

  BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB',
  END_CTAB_BLOCK = 'M  V30 END CTAB',

  BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM',
  END_ATOM_BLOCK = 'M  V30 END ATOM',

  BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND',
  END_BOND_BLOCK = 'M  V30 END BOND',
  BOND_CONFIG = ' CFG=',

  ATOM_TYPE_COL = 4,
  X_COL = 5,
  FIRST_BONDED_ATOM_COL = 5,
  BOND_TYPE_COL = 4,
  DUMMY_HEADER = `
     Datagrok

  0  0  0  0  0  0  0  0  0  0999 V3000`,
  BEGIN_CTAB = 'M  V30 BEGIN CTAB',
  COUNTS_LINE_START = 'M  V30 COUNTS ',
  COUNTS_LINE_DUMMY_END = ' 0 0 0',
  BEGIN_COLLECTION_BLOCK = 'M  V30 BEGIN COLLECTION',
  END_COLLECTION_BLOCK = 'M  V30 END COLLECTION',
  END_CTAB = 'M  V30 END CTAB',
  END = 'M  END',
}

export const EMPTY_MOL_V3000 = `
Empty input

  0  0  0  0  0  0            999 V3000
M  END`;
