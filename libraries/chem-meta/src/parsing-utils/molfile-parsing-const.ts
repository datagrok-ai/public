export enum MOLFILE_VERSION {
  V2000 = 'V2000',
  V3000 = 'V3000',
}

/** Constants relevant for parsing of Molfile V3K */
export const enum V3K {
  HEADER = 'V3000',
  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  /** Index shift from the beginning of the 'COUNTS' line to the number of atoms  */
  COUNTS_SHIFT = 14,
  BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM',
  BEGIN_BOND_BLOCK = 'M  V30 BEGIN BOND',
  END = 'M  END',
}

/** Constants relevant for parsing of Molfile V2K */
export enum V2K {
  HEADER = 'V2000',
  END = 'M  END',
}

