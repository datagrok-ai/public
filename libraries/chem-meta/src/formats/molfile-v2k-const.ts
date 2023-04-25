/** Constants related to the structure of Molfile V2K */
export enum V2K_CONST {
  HEADER = 'V2000',
  NUM_OF_HEADER_LINES = 3,

  NUM_OF_COUNTS_DIGITS = 3,
  ATOM_TYPE_COL = 4,
  FIRST_BONDED_ATOM_COL = 1,

  RGP_SHIFT = 8,
  RGP_LINE = 'M  RGP',
  A_LINE = 'A  ',

  END = 'M  END',
}

export const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;
