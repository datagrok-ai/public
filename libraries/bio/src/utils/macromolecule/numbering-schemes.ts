/** Static antibody numbering scheme region definitions.
 *
 * Contains FR/CDR boundaries for IMGT, Kabat, Chothia, and AHo schemes.
 * These are the canonical position ranges used after a numbering tool
 * (e.g. ANARCI) has assigned scheme-specific position names.
 */

export enum NumberingScheme {
  IMGT = 'IMGT',
  Kabat = 'Kabat',
  Chothia = 'Chothia',
  AHo = 'AHo',
}

export enum ChainType {
  Heavy = 'Heavy',
  Light_Kappa = 'Light_Kappa',
  Light_Lambda = 'Light_Lambda',
}

export interface SchemeRegionDef {
  name: string;
  type: 'FR' | 'CDR';
  startPosition: string;
  endPosition: string;
  chainType: ChainType;
}

/** IMGT region boundaries — the reference numbering scheme.
 * @see https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
 */
export const IMGT_REGIONS: Record<ChainType, SchemeRegionDef[]> = {
  [ChainType.Heavy]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '26', chainType: ChainType.Heavy},
    {name: 'CDR1', type: 'CDR', startPosition: '27', endPosition: '38', chainType: ChainType.Heavy},
    {name: 'FR2', type: 'FR', startPosition: '39', endPosition: '55', chainType: ChainType.Heavy},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '65', chainType: ChainType.Heavy},
    {name: 'FR3', type: 'FR', startPosition: '66', endPosition: '104', chainType: ChainType.Heavy},
    {name: 'CDR3', type: 'CDR', startPosition: '105', endPosition: '117', chainType: ChainType.Heavy},
    {name: 'FR4', type: 'FR', startPosition: '118', endPosition: '128', chainType: ChainType.Heavy},
  ],
  [ChainType.Light_Kappa]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '26', chainType: ChainType.Light_Kappa},
    {name: 'CDR1', type: 'CDR', startPosition: '27', endPosition: '38', chainType: ChainType.Light_Kappa},
    {name: 'FR2', type: 'FR', startPosition: '39', endPosition: '55', chainType: ChainType.Light_Kappa},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '65', chainType: ChainType.Light_Kappa},
    {name: 'FR3', type: 'FR', startPosition: '66', endPosition: '104', chainType: ChainType.Light_Kappa},
    {name: 'CDR3', type: 'CDR', startPosition: '105', endPosition: '117', chainType: ChainType.Light_Kappa},
    {name: 'FR4', type: 'FR', startPosition: '118', endPosition: '127', chainType: ChainType.Light_Kappa},
  ],
  [ChainType.Light_Lambda]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '26', chainType: ChainType.Light_Lambda},
    {name: 'CDR1', type: 'CDR', startPosition: '27', endPosition: '38', chainType: ChainType.Light_Lambda},
    {name: 'FR2', type: 'FR', startPosition: '39', endPosition: '55', chainType: ChainType.Light_Lambda},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '65', chainType: ChainType.Light_Lambda},
    {name: 'FR3', type: 'FR', startPosition: '66', endPosition: '104', chainType: ChainType.Light_Lambda},
    {name: 'CDR3', type: 'CDR', startPosition: '105', endPosition: '117', chainType: ChainType.Light_Lambda},
    {name: 'FR4', type: 'FR', startPosition: '118', endPosition: '127', chainType: ChainType.Light_Lambda},
  ],
};

/** Kabat region boundaries (differ mainly in CDR definitions). */
export const KABAT_REGIONS: Record<ChainType, SchemeRegionDef[]> = {
  [ChainType.Heavy]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '30', chainType: ChainType.Heavy},
    {name: 'CDR1', type: 'CDR', startPosition: '31', endPosition: '35', chainType: ChainType.Heavy},
    {name: 'FR2', type: 'FR', startPosition: '36', endPosition: '49', chainType: ChainType.Heavy},
    {name: 'CDR2', type: 'CDR', startPosition: '50', endPosition: '65', chainType: ChainType.Heavy},
    {name: 'FR3', type: 'FR', startPosition: '66', endPosition: '94', chainType: ChainType.Heavy},
    {name: 'CDR3', type: 'CDR', startPosition: '95', endPosition: '102', chainType: ChainType.Heavy},
    {name: 'FR4', type: 'FR', startPosition: '103', endPosition: '113', chainType: ChainType.Heavy},
  ],
  [ChainType.Light_Kappa]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '23', chainType: ChainType.Light_Kappa},
    {name: 'CDR1', type: 'CDR', startPosition: '24', endPosition: '34', chainType: ChainType.Light_Kappa},
    {name: 'FR2', type: 'FR', startPosition: '35', endPosition: '49', chainType: ChainType.Light_Kappa},
    {name: 'CDR2', type: 'CDR', startPosition: '50', endPosition: '56', chainType: ChainType.Light_Kappa},
    {name: 'FR3', type: 'FR', startPosition: '57', endPosition: '88', chainType: ChainType.Light_Kappa},
    {name: 'CDR3', type: 'CDR', startPosition: '89', endPosition: '97', chainType: ChainType.Light_Kappa},
    {name: 'FR4', type: 'FR', startPosition: '98', endPosition: '107', chainType: ChainType.Light_Kappa},
  ],
  [ChainType.Light_Lambda]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '23', chainType: ChainType.Light_Lambda},
    {name: 'CDR1', type: 'CDR', startPosition: '24', endPosition: '34', chainType: ChainType.Light_Lambda},
    {name: 'FR2', type: 'FR', startPosition: '35', endPosition: '49', chainType: ChainType.Light_Lambda},
    {name: 'CDR2', type: 'CDR', startPosition: '50', endPosition: '56', chainType: ChainType.Light_Lambda},
    {name: 'FR3', type: 'FR', startPosition: '57', endPosition: '88', chainType: ChainType.Light_Lambda},
    {name: 'CDR3', type: 'CDR', startPosition: '89', endPosition: '97', chainType: ChainType.Light_Lambda},
    {name: 'FR4', type: 'FR', startPosition: '98', endPosition: '107', chainType: ChainType.Light_Lambda},
  ],
};

/** Chothia region boundaries. */
export const CHOTHIA_REGIONS: Record<ChainType, SchemeRegionDef[]> = {
  [ChainType.Heavy]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '25', chainType: ChainType.Heavy},
    {name: 'CDR1', type: 'CDR', startPosition: '26', endPosition: '32', chainType: ChainType.Heavy},
    {name: 'FR2', type: 'FR', startPosition: '33', endPosition: '51', chainType: ChainType.Heavy},
    {name: 'CDR2', type: 'CDR', startPosition: '52', endPosition: '56', chainType: ChainType.Heavy},
    {name: 'FR3', type: 'FR', startPosition: '57', endPosition: '94', chainType: ChainType.Heavy},
    {name: 'CDR3', type: 'CDR', startPosition: '95', endPosition: '102', chainType: ChainType.Heavy},
    {name: 'FR4', type: 'FR', startPosition: '103', endPosition: '113', chainType: ChainType.Heavy},
  ],
  [ChainType.Light_Kappa]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '25', chainType: ChainType.Light_Kappa},
    {name: 'CDR1', type: 'CDR', startPosition: '26', endPosition: '32', chainType: ChainType.Light_Kappa},
    {name: 'FR2', type: 'FR', startPosition: '33', endPosition: '49', chainType: ChainType.Light_Kappa},
    {name: 'CDR2', type: 'CDR', startPosition: '50', endPosition: '52', chainType: ChainType.Light_Kappa},
    {name: 'FR3', type: 'FR', startPosition: '53', endPosition: '90', chainType: ChainType.Light_Kappa},
    {name: 'CDR3', type: 'CDR', startPosition: '91', endPosition: '96', chainType: ChainType.Light_Kappa},
    {name: 'FR4', type: 'FR', startPosition: '97', endPosition: '107', chainType: ChainType.Light_Kappa},
  ],
  [ChainType.Light_Lambda]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '25', chainType: ChainType.Light_Lambda},
    {name: 'CDR1', type: 'CDR', startPosition: '26', endPosition: '32', chainType: ChainType.Light_Lambda},
    {name: 'FR2', type: 'FR', startPosition: '33', endPosition: '49', chainType: ChainType.Light_Lambda},
    {name: 'CDR2', type: 'CDR', startPosition: '50', endPosition: '52', chainType: ChainType.Light_Lambda},
    {name: 'FR3', type: 'FR', startPosition: '53', endPosition: '90', chainType: ChainType.Light_Lambda},
    {name: 'CDR3', type: 'CDR', startPosition: '91', endPosition: '96', chainType: ChainType.Light_Lambda},
    {name: 'FR4', type: 'FR', startPosition: '97', endPosition: '107', chainType: ChainType.Light_Lambda},
  ],
};

/** AHo region boundaries. */
export const AHO_REGIONS: Record<ChainType, SchemeRegionDef[]> = {
  [ChainType.Heavy]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '24', chainType: ChainType.Heavy},
    {name: 'CDR1', type: 'CDR', startPosition: '25', endPosition: '40', chainType: ChainType.Heavy},
    {name: 'FR2', type: 'FR', startPosition: '41', endPosition: '55', chainType: ChainType.Heavy},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '78', chainType: ChainType.Heavy},
    {name: 'FR3', type: 'FR', startPosition: '79', endPosition: '108', chainType: ChainType.Heavy},
    {name: 'CDR3', type: 'CDR', startPosition: '109', endPosition: '138', chainType: ChainType.Heavy},
    {name: 'FR4', type: 'FR', startPosition: '139', endPosition: '149', chainType: ChainType.Heavy},
  ],
  [ChainType.Light_Kappa]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '24', chainType: ChainType.Light_Kappa},
    {name: 'CDR1', type: 'CDR', startPosition: '25', endPosition: '40', chainType: ChainType.Light_Kappa},
    {name: 'FR2', type: 'FR', startPosition: '41', endPosition: '55', chainType: ChainType.Light_Kappa},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '78', chainType: ChainType.Light_Kappa},
    {name: 'FR3', type: 'FR', startPosition: '79', endPosition: '108', chainType: ChainType.Light_Kappa},
    {name: 'CDR3', type: 'CDR', startPosition: '109', endPosition: '138', chainType: ChainType.Light_Kappa},
    {name: 'FR4', type: 'FR', startPosition: '139', endPosition: '149', chainType: ChainType.Light_Kappa},
  ],
  [ChainType.Light_Lambda]: [
    {name: 'FR1', type: 'FR', startPosition: '1', endPosition: '24', chainType: ChainType.Light_Lambda},
    {name: 'CDR1', type: 'CDR', startPosition: '25', endPosition: '40', chainType: ChainType.Light_Lambda},
    {name: 'FR2', type: 'FR', startPosition: '41', endPosition: '55', chainType: ChainType.Light_Lambda},
    {name: 'CDR2', type: 'CDR', startPosition: '56', endPosition: '78', chainType: ChainType.Light_Lambda},
    {name: 'FR3', type: 'FR', startPosition: '79', endPosition: '108', chainType: ChainType.Light_Lambda},
    {name: 'CDR3', type: 'CDR', startPosition: '109', endPosition: '138', chainType: ChainType.Light_Lambda},
    {name: 'FR4', type: 'FR', startPosition: '139', endPosition: '149', chainType: ChainType.Light_Lambda},
  ],
};

/** Lookup table: scheme → chain type → region defs */
export const SCHEME_REGIONS: Record<NumberingScheme, Record<ChainType, SchemeRegionDef[]>> = {
  [NumberingScheme.IMGT]: IMGT_REGIONS,
  [NumberingScheme.Kabat]: KABAT_REGIONS,
  [NumberingScheme.Chothia]: CHOTHIA_REGIONS,
  [NumberingScheme.AHo]: AHO_REGIONS,
};
