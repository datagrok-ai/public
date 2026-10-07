// Synced from Rusty-HMMER web/anarci/types.ts by web/sync-datagrok.mjs — do not edit here.
/** Shared types of the ANARCI port (ANARCI 2024.05.21, BSD-3-Clause,
 * Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, Claire
 * Marks; notice in licenses/ANARCI.txt). Python tuples become arrays with the
 * same element order, so values compare directly with the oracle's JSON. */

/** IMGT reference state (1-128) and its type: match, insert or delete. */
export type State = [number, 'm' | 'i' | 'd'];

/** One aligned column: the state and the 0-based sequence index (null at deletions). */
export type StateEntry = [State, number | null];

/** ANARCI's alignment of one domain to the 128 IMGT reference states. */
export type StateVector = StateEntry[];

/** A scheme position: number and insertion code (`' '` when none, else `A`, `B`, ... `ZZ`). */
export type Position = [number, string];

/** Numbered residues of one domain (`'-'` for an empty position). */
export type Numbering = [Position, string][];

/** `number_sequence_from_alignment` result: numbering, start and end (inclusive)
 * sequence indices (Python `None`, here null, when smoothing drops every state). */
export type NumberedDomain = [Numbering, number | null, number | null];

export type Scheme = 'imgt' | 'kabat' | 'chothia' | 'martin' | 'aho' | 'wolfguy';

/** ANARCI chain types: heavy, kappa, lambda, and TCR alpha, beta, gamma, delta. */
export type ChainType = 'H' | 'K' | 'L' | 'A' | 'B' | 'G' | 'D';

/** A Python `AssertionError` raised by ANARCI's numbering code. */
export class AnarciAssertion extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'AssertionError';
  }
}

/** Python `assert condition, message`. */
export function check(condition: unknown, message = ''): asserts condition {
  if (!condition) throw new AnarciAssertion(message);
}
