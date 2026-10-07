// Synced from Rusty-HMMER web/anarci/germline.ts by web/sync-datagrok.mjs — do not edit here.
/** ANARCI's germline assignment, ported from anarci.py (ANARCI 2024.05.21,
 * BSD-3-Clause, Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr
 * Kovaltsuk, Claire Marks; notice in licenses/ANARCI.txt): `get_identity` and
 * `run_germline_assignment`. Germline tables are loaded per chain type
 * (data/germlines-<chain>.json, generated from ANARCI's `all_germlines`). */

import type {Germlines} from './alignment.ts';
import type {StateVector} from './types.ts';

/** One chain type's germlines in ANARCI's dict order: `[species, [[gene, aligned 128-column sequence], ...]]`. */
export interface ChainGermlines {
  V: [string, [string, string][]][];
  J: [string, [string, string][]][];
}

/** `get_identity`: identity over the germline's non-gap columns. */
export function getIdentity(stateSequence: string, germline: string): number {
  if (!(stateSequence.length === germline.length && germline.length === 128))
    throw new Error('AssertionError: ');
  let n = 0;
  let m = 0;
  for (let i = 0; i < 128; i++) {
    if (germline[i] === '-') continue;
    if (stateSequence[i].toUpperCase() === germline[i]) m++;
    n++;
  }
  return n === 0 ? 0 : m / n;
}

/** Python `max(d, key=d.get)`: the first key with the largest value. */
function best(ids: [[string, string], number][]): [[string, string], number] {
  let top = ids[0];
  for (const entry of ids)
    if (entry[1] > top[1]) top = entry;
  return top;
}

/** `run_germline_assignment`. `table` is the chain's germlines or null when
 * ANARCI has none for the chain; `allSpecies` is ANARCI's `all_species`. */
export function runGermlineAssignment(stateVector: StateVector, sequence: string, table: ChainGermlines | null,
  allSpecies: string[], allowedSpecies: string[] | null): Germlines | Record<string, never> {
  const genes: Germlines = {v_gene: [null, null], j_gene: [null, null]};
  // dict(state_vector): the last sequence index stored for each match state.
  const matches = new Map<number, number | null>();
  for (const [[id, type], index] of stateVector)
    if (type === 'm') matches.set(id, index);
  let stateSequence = '';
  for (let i = 1; i <= 128; i++) {
    const index = matches.get(i);
    stateSequence += index === undefined || index === null ? '-' : sequence[index];
  }
  if (table === null || table.V.length === 0) return genes;
  const vBySpecies = new Map(table.V);
  let species = allowedSpecies;
  if (species !== null) {
    if (!species.every((s) => vBySpecies.has(s))) return {};
  } else
    species = allSpecies;

  const vIds: [[string, string], number][] = [];
  for (const s of species) {
    const genesOf = vBySpecies.get(s);
    if (!genesOf) continue;
    for (const [gene, germline] of genesOf)
      vIds.push([[s, gene], getIdentity(stateSequence, germline)]);
  }
  if (vIds.length === 0) throw new Error('ValueError: max() arg is an empty sequence');
  genes.v_gene = best(vIds);
  const vSpecies = genes.v_gene[0][0];
  const jGenes = new Map(table.J).get(vSpecies);
  if (table.J.length > 0 && jGenes) {
    const jIds: [[string, string], number][] = jGenes.map(([gene, germline]) =>
      [[vSpecies, gene], getIdentity(stateSequence, germline)]);
    if (jIds.length === 0) throw new Error('ValueError: max() arg is an empty sequence');
    genes.j_gene = best(jIds);
  }
  return genes;
}
