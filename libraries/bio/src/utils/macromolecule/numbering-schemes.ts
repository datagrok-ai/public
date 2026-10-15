/** Static antibody numbering scheme region definitions.
 *
 * Contains FR/CDR boundaries for IMGT, Kabat, Chothia, Martin, AHo and Wolfguy
 * schemes (immunoglobulin chains) and IMGT/AHo for T-cell receptor chains.
 * These are the canonical position ranges used after a numbering tool
 * (e.g. ANARCI) has assigned scheme-specific position names.
 */

export enum NumberingScheme {
  IMGT = 'IMGT',
  Kabat = 'Kabat',
  Chothia = 'Chothia',
  Martin = 'Martin',
  AHo = 'AHo',
  Wolfguy = 'Wolfguy',
}

export enum ChainType {
  Heavy = 'Heavy',
  Light_Kappa = 'Light_Kappa',
  Light_Lambda = 'Light_Lambda',
  TCR_Alpha = 'TCR_Alpha',
  TCR_Beta = 'TCR_Beta',
  TCR_Gamma = 'TCR_Gamma',
  TCR_Delta = 'TCR_Delta',
}

/** Region definitions of each chain type a scheme numbers (immunoglobulin
 *  only for Kabat, Chothia, Martin and Wolfguy). */
export type SchemeRegions = Partial<Record<ChainType, SchemeRegionDef[]>>;

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
export const IMGT_REGIONS: SchemeRegions = {
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
export const KABAT_REGIONS: SchemeRegions = {
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
export const CHOTHIA_REGIONS: SchemeRegions = {
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

/** FR1, CDR1, FR2, CDR2, FR3, CDR3 and FR4 of `chainType` from consecutive bounds. */
function regionList(chainType: ChainType, bounds: [number, number][]): SchemeRegionDef[] {
  const names = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'];
  return bounds.map(([start, end], i) => ({name: names[i], type: i % 2 === 1 ? 'CDR' : 'FR',
    startPosition: String(start), endPosition: String(end), chainType}));
}

/** IMGT numbers all V domains alike (the IMGT unique numbering), T-cell
 *  receptor chains included; VJ (alpha, gamma) and VDJ (beta, delta) chains
 *  share the immunoglobulin boundaries, FR4 up to 128. */
const IMGT_TCR: [number, number][] = [[1, 26], [27, 38], [39, 55], [56, 65], [66, 104], [105, 117], [118, 128]];
for (const chain of [ChainType.TCR_Alpha, ChainType.TCR_Beta, ChainType.TCR_Gamma, ChainType.TCR_Delta])
  IMGT_REGIONS[chain] = regionList(chain, IMGT_TCR);

/** Martin (enhanced Chothia) numbering differs from Chothia in framework
 *  insertion positions only; its CDRs are the Chothia ones. */
export const MARTIN_REGIONS: SchemeRegions = Object.fromEntries(Object.entries(CHOTHIA_REGIONS).map(
  ([chain, defs]) => [chain, defs!.map((d) => ({...d}))]));

/** Wolfguy encodes the region in the position number: hundreds per region
 *  (heavy 1xx-4xx, light 5xx-8xx), framework positions x01-x50, CDR positions x51-x99. */
function wolfguyRegions(chainType: ChainType, first: number): SchemeRegionDef[] {
  return regionList(chainType, Array.from({length: 7}, (_, i): [number, number] => {
    const hundred = (first + Math.floor(i / 2)) * 100;
    return i % 2 === 1 ? [hundred + 51, hundred + 99] : [hundred + 1, hundred + 50];
  }));
}

export const WOLFGUY_REGIONS: SchemeRegions = {
  [ChainType.Heavy]: wolfguyRegions(ChainType.Heavy, 1),
  [ChainType.Light_Kappa]: wolfguyRegions(ChainType.Light_Kappa, 5),
  [ChainType.Light_Lambda]: wolfguyRegions(ChainType.Light_Lambda, 5),
};

/** AHo region boundaries. */
export const AHO_REGIONS: SchemeRegions = {
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
/** AHo numbers T-cell receptor V domains on the same 149 positions. */
const AHO_TCR: [number, number][] = [[1, 24], [25, 40], [41, 55], [56, 78], [79, 108], [109, 138], [139, 149]];
for (const chain of [ChainType.TCR_Alpha, ChainType.TCR_Beta, ChainType.TCR_Gamma, ChainType.TCR_Delta])
  AHO_REGIONS[chain] = regionList(chain, AHO_TCR);

export const SCHEME_REGIONS: Record<NumberingScheme, SchemeRegions> = {
  [NumberingScheme.IMGT]: IMGT_REGIONS,
  [NumberingScheme.Kabat]: KABAT_REGIONS,
  [NumberingScheme.Chothia]: CHOTHIA_REGIONS,
  [NumberingScheme.Martin]: MARTIN_REGIONS,
  [NumberingScheme.AHo]: AHO_REGIONS,
  [NumberingScheme.Wolfguy]: WOLFGUY_REGIONS,
};
