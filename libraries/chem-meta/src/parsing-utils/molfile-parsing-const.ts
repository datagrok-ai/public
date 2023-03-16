export enum MOLFILE_VERSION {
  V2000 = 'V2000',
  V3000 = 'V3000',
}

/** Constants relevant for parsing of Molfile V3K */
export const enum V3K_CONST {
  HEADER = 'V3000',
  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  /** Index shift from the beginning of the 'COUNTS' line to the number of atoms  */
  COUNTS_SHIFT = 14,
  BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM',
  BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND',
  END = 'M  END',
  ATOM_TYPE_COL = 4,
  X_COL = 5,
  FIRST_BONDED_ATOM_COL = 5,
}

/** Constants relevant for parsing of Molfile V2K */
export enum V2K_CONST {
  HEADER = 'V2000',
  END = 'M  END',
  NUM_OF_HEADER_LINES = 3,
  NUM_OF_COUNTS_DIGITS = 3,
  ATOM_TYPE_COL = 4,
  FIRST_BONDED_ATOM_COL = 1,
}

