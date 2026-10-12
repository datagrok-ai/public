// Synced from Rusty-HMMER web/anarci/anarci.ts by web/sync-datagrok.mjs — do not edit here.
/** ANARCI (Antibody Numbering and Antigen Receptor ClassIfication) on the
 * HMMER engine. Port of anarci.py `anarci()` and
 * `number_sequences_from_alignment` (ANARCI 2024.05.21, BSD-3-Clause,
 * Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, Claire
 * Marks; notice in licenses/ANARCI.txt), verified against unmodified ANARCI
 * (oracle/anarci).
 *
 * One deliberate difference: ANARCI raises (and loses the whole batch) when
 * numbering one sequence fails; here that sequence gets an `error` and the
 * others are numbered. */

import type {HmmDatabase} from '../engine.ts';
import {checkForJ, type DomainDetails, type HitRow, runHmmer, type SequenceAlignment} from './alignment.ts';
import {type ChainGermlines, runGermlineAssignment} from './germline.ts';
import {numberSequenceFromAlignment, validateNumbering} from './schemes.ts';
import type {ChainType, NumberedDomain, Scheme} from './types.ts';

export type {DomainDetails, HitRow, SequenceAlignment} from './alignment.ts';
export type {ChainGermlines} from './germline.ts';

export const ALL_CHAINS: ChainType[] = ['H', 'K', 'L', 'A', 'B', 'G', 'D'];
export const IG_CHAINS: ChainType[] = ['H', 'K', 'L'];
/** ANARCI's default `allowed_species` (`anarci()`, `number()` and the command line). */
export const DEFAULT_SPECIES = ['human', 'mouse'];
export const SCHEMES: Scheme[] = ['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy'];

/** ANARCI's model database and data tables. */
export interface AnarciData {
  database: HmmDatabase;
  /** `get_hmm_length` per `species_chain` (data/hmm-lengths.json). */
  hmmLengths: Record<string, number>;
  /** `all_species` (data/species.json). */
  allSpecies: string[];
  /** Germline tables per chain type; only needed with `assignGermline`. */
  germlines?: (chain: string) => ChainGermlines | null;
}

export interface AnarciOptions {
  scheme: Scheme;
  /** Chain types to number. Default: all for IMGT and AHo, H/K/L otherwise
   * (as the ANARCI command line does; ANARCI cannot number TCRs in other schemes). */
  allow?: Iterable<string>;
  /** `allowed_species`: preferred species for the HMM hit and the germline
   * assignment; `null` allows all. Default `['human', 'mouse']`. */
  allowedSpecies?: string[] | null;
  assignGermline?: boolean;
  bitScoreThreshold?: number;
}

export interface SequenceResult {
  /** Numbered domains, or null when none was numbered (ANARCI's `numbered[i]`). */
  numbered: NumberedDomain[] | null;
  /** Details of the numbered domains (ANARCI's `alignment_details[i]`). */
  details: DomainDetails[] | null;
  /** All hits over the threshold (ANARCI's `hit_tables[i]`, without its header row). */
  hitTable: HitRow[];
  /** Why this sequence could not be numbered (ANARCI would raise). */
  error?: string;
}

/** `number_sequences_from_alignment` for one sequence. Clones the alignment's
 * details, so one alignment can be numbered with several schemes. */
export function numberFromAlignment(data: AnarciData, name: string, sequence: string, alignment: SequenceAlignment,
  options: AnarciOptions): SequenceResult {
  const scheme = options.scheme;
  const allow = new Set(options.allow ?? (scheme === 'imgt' || scheme === 'aho' ? ALL_CHAINS : IG_CHAINS));
  const species = options.allowedSpecies === undefined ? DEFAULT_SPECIES : options.allowedSpecies;
  const numbered: NumberedDomain[] = [];
  const details: DomainDetails[] = [];
  try {
    for (let di = 0; di < alignment.stateVectors.length; di++) {
      const vector = alignment.stateVectors[di];
      const detail: DomainDetails = {...alignment.details[di], scheme, query_name: name};
      if (vector.length > 0 && allow.has(detail.chain_type)) {
        numbered.push(validateNumbering(
          numberSequenceFromAlignment(vector, sequence, scheme, detail.chain_type), name, sequence));
        if (options.assignGermline) {
          const table = data.germlines ? data.germlines(detail.chain_type) : null;
          detail.germlines = runGermlineAssignment(vector, sequence, table, data.allSpecies, species);
        }
        details.push(detail);
      }
    }
  } catch (e) {
    const error = e instanceof Error ? `${e.name}: ${e.message}` : String(e);
    return {numbered: null, details: null, hitTable: alignment.hitTable, error};
  }
  return numbered.length > 0 ?
    {numbered, details, hitTable: alignment.hitTable} :
    {numbered: null, details: null, hitTable: alignment.hitTable};
}

/** `run_hmmer` followed by `check_for_j`: the scheme-independent alignments. */
export function align(data: AnarciData, sequences: [string, string][], options: Omit<AnarciOptions, 'scheme'> = {}):
  SequenceAlignment[] {
  const species = options.allowedSpecies === undefined ? DEFAULT_SPECIES : options.allowedSpecies;
  const alignments = runHmmer(data.database, sequences, data.hmmLengths, options.bitScoreThreshold ?? 80, species);
  checkForJ(data.database, sequences, alignments, data.hmmLengths);
  return alignments;
}

/** `anarci()`: identify, number and (optionally) germline-assign the domains of
 * `(name, sequence)` pairs. */
export function anarci(data: AnarciData, sequences: [string, string][], options: AnarciOptions): SequenceResult[] {
  const alignments = align(data, sequences, options);
  return sequences.map(([name, sequence], i) => numberFromAlignment(data, name, sequence, alignments[i], options));
}
