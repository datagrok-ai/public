// Synced from Rusty-HMMER web/anarci/alignment.ts by web/sync-datagrok.mjs — do not edit here.
/** ANARCI's alignment step, ported from anarci.py (ANARCI 2024.05.21,
 * BSD-3-Clause, Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr
 * Kovaltsuk, Claire Marks; notice in licenses/ANARCI.txt): `run_hmmer`,
 * `_parse_hmmer_query`, `_domains_are_same`, `_hmm_alignment_to_states` and
 * `check_for_j`.
 *
 * ANARCI runs C `hmmscan` and parses its text report with Biopython; here the
 * same comparisons run in the HMMER engine, and each reported domain becomes
 * the HSP Biopython would build: coordinates from the domain table (0-based
 * starts), scores and E-values as printed (`%6.1f`, `%5.1f`, `%9.2g`) and read
 * back, and the alignment's RF and PP lines. */

import type {HmmDatabase} from '../engine.ts';
import type {StateEntry, StateVector} from './types.ts';

/** One `hit_table` row: id, description, evalue, bitscore, bias, query_start, query_end. */
export type HitRow = [string, string, number, number, number, number, number];

/** The per-domain description ANARCI builds (`top_descriptions`), later extended
 * by numbering with `scheme`, `query_name` and `germlines`. */
export interface DomainDetails {
  id: string;
  description: string;
  evalue: number;
  bitscore: number;
  bias: number;
  query_start: number | null;
  query_end: number;
  species: string;
  chain_type: string;
  scheme?: string;
  query_name?: string;
  germlines?: Germlines | Record<string, never>;
}

/** `run_germline_assignment` result: `[[species, gene], identity]` or `[null, null]`. */
export interface Germlines {
  v_gene: [[string, string], number] | [null, null];
  j_gene: [[string, string], number] | [null, null];
}

/** One sequence's `run_hmmer` result: (hit_table without its header row, state vectors, details). */
export interface SequenceAlignment {
  hitTable: HitRow[];
  stateVectors: StateVector[];
  details: DomainDetails[];
}

/** The Biopython HSP fields ANARCI reads. */
interface Hsp {
  hitId: string;
  hitDescription: string;
  evalue: number;
  bitscore: number;
  bias: number;
  queryStart: number;
  queryEnd: number;
  hitStart: number;
  hitEnd: number;
  reference: string;
  posterior: string;
  order?: number;
}

/** `all_reference_states`: the 128 IMGT match states. */
const REFERENCE_STATES = Array.from({length: 128}, (_, i) => i + 1);

/** Python `list[index]` for a non-negative index, raising IndexError. */
function at<T>(list: T[], index: number): T {
  if (index < 0 || index >= list.length) throw new Error('IndexError: list index out of range');
  return list[index];
}

/** `_domains_are_same`: whether two HSPs overlap on the query. */
function domainsAreSame(a: Hsp, b: Hsp): boolean {
  const [first, second] = a.queryStart <= b.queryStart ? [a, b] : [b, a];
  return !(second.queryStart >= first.queryEnd);
}

/** `get_hmm_length(species, ctype)`. */
function hmmLength(lengths: Record<string, number>, id: string): number {
  return lengths[id] ?? 128;
}

/** `_hmm_alignment_to_states`. */
function alignmentToStates(hsp: Hsp, n: number, seqLength: number, lengths: Record<string, number>): StateVector {
  let reference = hsp.reference;
  let states = hsp.posterior;
  if (reference.length !== states.length) {
    throw new Error('AssertionError: Aligned reference and state strings had different lengths. ' +
      'Don\'t know how to handle');
  }
  let hmmStart = hsp.hitStart;
  let hmmEnd = hsp.hitEnd;
  let seqStart = hsp.queryStart;
  let seqEnd = hsp.queryEnd;
  const parts = hsp.hitId.split('_');
  if (parts.length !== 2) throw new Error('ValueError: too many values to unpack');
  const fullLength = hmmLength(lengths, hsp.hitId);

  // Up to 5 unmatched N-terminal states of the first domain are numbered.
  if (hsp.order === 0 && hmmStart && hmmStart < 5) {
    let extend = hmmStart;
    if (hmmStart > seqStart)
      extend = Math.min(seqStart, hmmStart - seqStart);
    states = '8'.repeat(Math.max(0, extend)) + states;
    reference = 'x'.repeat(Math.max(0, extend)) + reference;
    seqStart -= extend;
    hmmStart -= extend;
  }
  // Extend a single domain to the end of the J region.
  if (n === 1 && seqEnd < seqLength && (123 < hmmEnd && hmmEnd < fullLength)) {
    const extend = Math.min(fullLength - hmmEnd, seqLength - seqEnd);
    states = states + '8'.repeat(Math.max(0, extend));
    reference = reference + 'x'.repeat(Math.max(0, extend));
    seqEnd += extend;
    hmmEnd += extend;
  }

  const hmmStates = REFERENCE_STATES.slice(Math.max(0, hmmStart), Math.max(0, hmmEnd));
  const sequenceIndices: number[] = [];
  for (let i = seqStart; i < seqEnd; i++) sequenceIndices.push(i);
  let h = 0;
  let s = 0;
  const vector: StateVector = [];
  for (let i = 0; i < states.length; i++) {
    let type: 'm' | 'i' | 'd' = reference[i] === 'x' ? 'm' : 'i';
    let index: number | null;
    if (states[i] === '.') {
      type = 'd';
      index = null;
    } else
      index = at(sequenceIndices, s);

    vector.push([[at(hmmStates, h), type], index]);
    if (type === 'm') {
      h++;
      s++;
    } else if (type === 'i')
      s++;
    else
      h++;
  }
  return vector;
}

/** The HSPs Biopython yields for one query: reported domains of reported hits, in report order. */
function hsps(database: HmmDatabase, result: ReturnType<HmmDatabase['scan']>[number]): Hsp[] {
  const out: Hsp[] = [];
  for (const hit of result.hits) {
    const model = database.info.models[hit.target];
    for (const domain of hit.domains) {
      if (!domain.reported || !domain.alignment) continue;
      out.push({
        hitId: model.name,
        hitDescription: model.description,
        evalue: domain.printed.iEvalue,
        bitscore: domain.printed.bitscore,
        bias: domain.printed.bias,
        queryStart: domain.aliFrom - 1,
        queryEnd: domain.aliTo,
        hitStart: domain.hmmFrom - 1,
        hitEnd: domain.hmmTo,
        reference: domain.alignment.reference ?? '',
        posterior: domain.alignment.posterior ?? '',
      });
    }
  }
  return out;
}

/** `_parse_hmmer_query`. */
function parseQuery(found: Hsp[], seqLength: number, threshold: number, hmmerSpecies: string[] | null,
  lengths: Record<string, number>): SequenceAlignment {
  const hitTable: HitRow[] = [];
  let domains: Hsp[] = [];
  let descriptions: DomainDetails[] = [];
  const stateVectors: StateVector[] = [];
  if (found.length > 0) {
    let list = found;
    if (hmmerSpecies && hmmerSpecies.length > 0) {
      const correct: Hsp[] = [];
      for (const hsp of found) {
        if (hsp.bitscore >= threshold) {
          for (const species of hmmerSpecies)
            if (hsp.hitId.startsWith(species)) correct.push(hsp);
        }
      }
      list = correct.length > 0 ? correct : found;
    }
    // Array.prototype.sort is stable, as Python's sorted().
    for (const hsp of [...list].sort((a, b) => a.evalue - b.evalue)) {
      if (hsp.bitscore >= threshold) {
        let isNew = true;
        for (const domain of domains) {
          if (domainsAreSame(domain, hsp)) {
            isNew = false;
            break;
          }
        }
        const row: HitRow = [hsp.hitId, hsp.hitDescription, hsp.evalue, hsp.bitscore, hsp.bias, hsp.queryStart,
          hsp.queryEnd];
        hitTable.push(row);
        if (isNew) {
          domains.push(hsp);
          descriptions.push({
            id: row[0], description: row[1], evalue: row[2], bitscore: row[3], bias: row[4],
            query_start: row[5], query_end: row[6], species: '', chain_type: '',
          });
        }
      }
    }
    const ordering = domains.map((_, i) => i).sort((a, b) => domains[a].queryStart - domains[b].queryStart);
    domains = ordering.map((i) => domains[i]);
    descriptions = ordering.map((i) => descriptions[i]);
  }
  for (let i = 0; i < domains.length; i++) {
    domains[i].order = i;
    const [species, chain] = descriptions[i].id.split('_');
    const vector = alignmentToStates(domains[i], domains.length, seqLength, lengths);
    stateVectors.push(vector);
    // Python: dict keys keep their first position; species and chain_type are appended.
    descriptions[i].species = species;
    descriptions[i].chain_type = chain;
    descriptions[i].query_start = vector[0][1];
  }
  return {hitTable, stateVectors, details: descriptions};
}

/** `run_hmmer`: scan every sequence against the germline database and parse
 * each query. `hmmerSpecies` prefers hits of those species when any reaches
 * the threshold (ANARCI's `allowed_species`). */
export function runHmmer(database: HmmDatabase, sequences: [string, string][], lengths: Record<string, number>,
  threshold = 80, hmmerSpecies: string[] | null = null): SequenceAlignment[] {
  const results = database.scan(sequences.map(([name, residues]) => ({name, residues})), {alignments: true});
  return results.map((result, i) => {
    if (result.status !== 0 && result.status !== 2)
      throw new Error(`hmmscan failed for ${sequences[i][0]} (status ${result.status})`);
    return parseQuery(hsps(database, result), sequences[i][1].length, threshold, hmmerSpecies, lengths);
  });
}

/** Python `dict(pairs).get(key)`: the last value stored under the key. */
function lookupState(vector: StateVector, id: number, type: string): number | null | undefined {
  let found: number | null | undefined;
  for (const [[stateId, stateType], index] of vector)
    if (stateId === id && stateType === type) found = index;
  return found;
}

/** `check_for_j`: when a single-domain alignment stops before the J region and
 * much sequence remains (a long CDR3), look for the J region after the
 * conserved cysteine (IMGT 104) with a bit score threshold of 10 and splice
 * the CDR3 between the V and J regions. Modifies `alignments` in place. */
export function checkForJ(database: HmmDatabase, sequences: [string, string][], alignments: SequenceAlignment[],
  lengths: Record<string, number>): void {
  for (let i = 0; i < sequences.length; i++) {
    const alignment = alignments[i];
    if (alignment.stateVectors.length !== 1) continue;
    const ali = alignment.stateVectors[0];
    const last = ali[ali.length - 1];
    const lastState = last[0][0];
    const lastIndex = last[1];
    if (lastState >= 120) continue;
    if (lastIndex === null) throw new Error('TypeError: unsupported operand type(s) for +: \'NoneType\' and \'int\'');
    if (!(lastIndex + 30 < sequences[i][1].length)) continue;
    const cys = lookupState(ali, 104, 'm');
    if (cys === undefined || cys === null) continue;
    const cysAt = ali.findIndex(([[id, type], index]) => id === 104 && type === 'm' && index === cys);
    const [rescan] = runHmmer(database, [[sequences[i][0], sequences[i][1].slice(cys + 1)]], lengths, 10);
    const states = rescan.stateVectors;
    if (!(states.length > 0 && states[0][states[0].length - 1][0][0] >= 126 && states[0][0][0][0] <= 117)) continue;
    const vRegion = ali.slice(0, cysAt + 1);
    const jRegion: StateEntry[] = [];
    for (const [state, index] of states[0])
      if (state[0] >= 117 && index !== null) jRegion.push([state, index + cys + 1]);
    const first = at(jRegion, 0)[1] as number;
    const cdrRegion: StateEntry[] = [];
    let next = 105;
    for (let si = cys + 1; si < first; si++) {
      if (next >= 116)
        cdrRegion.push([[116, 'i'], si]);
      else {
        cdrRegion.push([[next, 'm'], si]);
        next++;
      }
    }
    alignment.stateVectors[0] = [...vRegion, ...cdrRegion, ...jRegion];
    alignment.details[0].query_end = (jRegion[jRegion.length - 1][1] as number) + 1;
  }
}
