/** Constants related to the structure of Molfile V3K */
export const enum V3K_CONST {
  HEADER = 'V3000',
  HEADER_SECOND_LINE = '  0  0  0  0  0  0            999 V3000\n',

  BEGIN_DATA_LINE = 'M  V30 ',

  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  /** Index shift from the beginning of the 'COUNTS' line to the number of atoms  */
  COUNTS_SHIFT = 14,

  BEGIN_CTAB_BLOCK = 'M  V30 BEGIN CTAB',
  END_CTAB_BLOCK = 'M  V30 END CTAB',

  BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM',
  // /** Shift from the begginning of bond/atom block line to the corresponding
  //  * index */
  // IDX_SHIFT = 7,
  END_ATOM_BLOCK = 'M  V30 END ATOM',

  BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND',
  END_BOND_BLOCK = 'M  V30 END BOND',
  BOND_CONFIG = ' CFG=',

  ATOM_TYPE_COL = 4,
  X_COL = 5,
  FIRST_BONDED_ATOM_COL = 5,

  END = 'M  END',
}

export const EMPTY_MOL_V3000 = `
Empty input

  0  0  0  0  0  0            999 V3000
M  END`;
