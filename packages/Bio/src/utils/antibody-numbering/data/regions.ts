/**
 * Region definitions per scheme, keyed by chain group.
 * Each region: [name, startPosition, endPosition].
 * Mirrors the Python number_antibody.py definitions.
 */
import type { Scheme, ChainGroup } from '../types';

export type RegionDef = [name: string, start: number, end: number];

const IMGT_REGIONS: Record<ChainGroup, RegionDef[]> = {
  Heavy: [
    ['FR1', 1, 26], ['CDR1', 27, 38], ['FR2', 39, 55],
    ['CDR2', 56, 65], ['FR3', 66, 104], ['CDR3', 105, 117], ['FR4', 118, 128],
  ],
  Light: [
    ['FR1', 1, 26], ['CDR1', 27, 38], ['FR2', 39, 55],
    ['CDR2', 56, 65], ['FR3', 66, 104], ['CDR3', 105, 117], ['FR4', 118, 127],
  ],
};

const KABAT_REGIONS: Record<ChainGroup, RegionDef[]> = {
  Heavy: [
    ['FR1', 1, 30], ['CDR1', 31, 35], ['FR2', 36, 49],
    ['CDR2', 50, 65], ['FR3', 66, 94], ['CDR3', 95, 102], ['FR4', 103, 113],
  ],
  Light: [
    ['FR1', 1, 23], ['CDR1', 24, 34], ['FR2', 35, 49],
    ['CDR2', 50, 56], ['FR3', 57, 88], ['CDR3', 89, 97], ['FR4', 98, 107],
  ],
};

const CHOTHIA_REGIONS: Record<ChainGroup, RegionDef[]> = {
  Heavy: [
    ['FR1', 1, 25], ['CDR1', 26, 32], ['FR2', 33, 51],
    ['CDR2', 52, 56], ['FR3', 57, 94], ['CDR3', 95, 102], ['FR4', 103, 113],
  ],
  Light: [
    ['FR1', 1, 25], ['CDR1', 26, 32], ['FR2', 33, 49],
    ['CDR2', 50, 52], ['FR3', 53, 90], ['CDR3', 91, 96], ['FR4', 97, 107],
  ],
};

const AHO_REGIONS: Record<ChainGroup, RegionDef[]> = {
  Heavy: [
    ['FR1', 1, 24], ['CDR1', 25, 40], ['FR2', 41, 55],
    ['CDR2', 56, 78], ['FR3', 79, 108], ['CDR3', 109, 138], ['FR4', 139, 149],
  ],
  Light: [
    ['FR1', 1, 24], ['CDR1', 25, 40], ['FR2', 41, 55],
    ['CDR2', 56, 78], ['FR3', 79, 108], ['CDR3', 109, 138], ['FR4', 139, 149],
  ],
};

const SCHEME_REGIONS: Record<Scheme, Record<ChainGroup, RegionDef[]>> = {
  imgt: IMGT_REGIONS,
  kabat: KABAT_REGIONS,
  chothia: CHOTHIA_REGIONS,
  aho: AHO_REGIONS,
};

export function getRegions(scheme: Scheme, chainGroup: ChainGroup): RegionDef[] {
  return SCHEME_REGIONS[scheme]?.[chainGroup] ?? IMGT_REGIONS.Heavy;
}
