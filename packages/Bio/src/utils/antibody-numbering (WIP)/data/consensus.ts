/**
 * Central consensus profile registry.
 * Maps (scheme, chainType) -> profile data.
 */
import type { Scheme, ChainType } from '../types';
import { IMGT_H, IMGT_K, IMGT_L, CTERM_H, CTERM_K, CTERM_L } from './consensus-imgt';
import { KABAT_H, KABAT_K, KABAT_L } from './consensus-kabat';
import { MARTIN_H, MARTIN_K, MARTIN_L } from './consensus-martin';
import { AHO_H, AHO_K, AHO_L } from './consensus-aho';

export type ProfileData = [number, string[]][];

/** Get the consensus profile for a given scheme and chain type */
export function getConsensusProfile(scheme: Scheme, chain: ChainType): ProfileData {
  const key = `${scheme}_${chain}`;
  const profiles: Record<string, ProfileData> = {
    imgt_H: IMGT_H, imgt_K: IMGT_K, imgt_L: IMGT_L,
    kabat_H: KABAT_H, kabat_K: KABAT_K, kabat_L: KABAT_L,
    chothia_H: MARTIN_H, chothia_K: MARTIN_K, chothia_L: MARTIN_L,
    aho_H: AHO_H, aho_K: AHO_K, aho_L: AHO_L,
  };
  return profiles[key] ?? IMGT_H;
}

/** Get the C-terminal finder profile for a chain type */
export function getCtermProfile(chain: ChainType): ProfileData {
  const profiles: Record<ChainType, ProfileData> = {
    H: CTERM_H, K: CTERM_K, L: CTERM_L,
  };
  return profiles[chain];
}

export { IMGT_H, IMGT_K, IMGT_L, CTERM_H, CTERM_K, CTERM_L };
export { KABAT_H, KABAT_K, KABAT_L };
export { MARTIN_H, MARTIN_K, MARTIN_L };
export { AHO_H, AHO_K, AHO_L };
