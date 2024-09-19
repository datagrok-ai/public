import {PolymerTypes} from '../helm/consts';

export const monomerWorksConsts = {
// constants for parsing molfile V2000
  V2K_RGP_SHIFT: 8,
  V2K_RGP_LINE: 'M  RGP',
  V2K_A_LINE: 'A  ',
  // constants for parsing/reconstruction of molfile V3000
  V3K_COUNTS_SHIFT: 14,
  V3K_IDX_SHIFT: 7,
  V3K_HEADER_FIRST_LINE: '\nDatagrok macromolecule handler\n\n',
  V3K_HEADER_SECOND_LINE: '  0  0  0  0  0  0            999 V3000\n',
  V3K_BEGIN_CTAB_BLOCK: 'M  V30 BEGIN CTAB\n',
  V3K_END_CTAB_BLOCK: 'M  V30 END CTAB\n',
  V3K_BEGIN_COUNTS_LINE: 'M  V30 COUNTS ',
  V3K_COUNTS_LINE_ENDING: ' 0 0 0\n',
  V3K_BEGIN_ATOM_BLOCK: 'M  V30 BEGIN ATOM\n',
  V3K_END_ATOM_BLOCK: 'M  V30 END ATOM\n',
  V3K_BEGIN_BOND_BLOCK: 'M  V30 BEGIN BOND\n',
  V3K_END_BOND_BLOCK: 'M  V30 END BOND\n',
  V3K_BOND_CONFIG: ' CFG=',
  V3K_BEGIN_DATA_LINE: 'M  V30 ',
  V3K_END: 'M  END',
  PRECISION_FACTOR: 10_000, // HELMCoreLibrary has 4 significant digits after decimal point in atom coordinates
  // symbols for the corresponding monomers in HELM library

  DEOXYRIBOSE: {polymerType: PolymerTypes.RNA, symbol: 'd'},
  RIBOSE: {polymerType: PolymerTypes.RNA, symbol: 'r'},
  PHOSPHATE: {polymerType: PolymerTypes.RNA, symbol: 'p'},

  OXYGEN: 'O',
  HYDROGEN: 'H',
} as const;

