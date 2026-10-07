// Synced from Rusty-HMMER web/search.ts by web/sync-datagrok.mjs — do not edit here.
/** `hmmsearch` split across workers.
 *
 * Every target sequence is compared independently, so a search can run in
 * parts: each part gets the total number of sequences as Z (`-Z`) and returns
 * every registered hit (`allHits`). The parts are then merged and thresholded
 * here exactly as `p7_tophits_Threshold` (HMMER 3.4 `src/p7_tophits.c`,
 * BSD-3-Clause, notice in licenses/HMMER.txt) does for one search: targets by
 * E-value (P-value × Z) or score, domZ as the number of reported targets,
 * domains by conditional E-value (P-value × domZ) or score, and the bug #h74
 * workaround for identical domain alignments. P-values are the engine's own
 * `exp(lnP)`, so every product is bit-identical to the unsplit search. */

import {HitFlags, type Hit, type QueryResult, type SearchOptions} from './engine.ts';

/** Options for each part of a search over `total` sequences. */
export function partOptions(options: SearchOptions, total: number): SearchOptions {
  return {...options, allHits: true, z: options.z ?? total};
}

/** A part's result with the index of its first sequence in the whole search. */
export interface SearchPart {
  result: QueryResult;
  offset: number;
}

export interface MergedSearch {
  /** All registered hits, ranked; target indices refer to the whole search. */
  hits: Hit[];
  z: number;
  domZ: number;
  reported: number;
  included: number;
}

/** Merge parts run with {@link partOptions}. `names` (one per sequence) break
 * ranking ties as C's `strcmp` does; by default sequences are named by index. */
export function mergeSearch(parts: SearchPart[], options: SearchOptions, total: number,
  names?: string[]): MergedSearch {
  const z = options.z ?? total;
  const hits: Hit[] = parts.flatMap(({result, offset}) =>
    result.hits.map((hit) => ({...hit, target: hit.target + offset, domains: hit.domains.map((d) => ({...d}))})));
  const cutoffs = options.cutoff !== undefined;
  const targetE = options.evalue ?? 10;
  const incE = options.incEvalue ?? 0.01;
  const domE = options.domEvalue ?? 10;
  const incdomE = options.incdomEvalue ?? 0.01;
  const both = HitFlags.reported | HitFlags.included;

  // Without model cutoffs, target flags are set here; with them, the pipeline set them.
  if (!cutoffs) {
    for (const hit of hits) {
      hit.flags &= ~both;
      const reportable = options.score === undefined ? hit.p * z <= targetE : hit.score >= options.score;
      if (!reportable) continue;
      hit.flags |= HitFlags.reported;
      const includable = options.incScore === undefined ? hit.p * z <= incE : hit.score >= options.incScore;
      if (includable) hit.flags |= HitFlags.included;
    }
  }
  const reported = hits.filter((h) => h.flags & HitFlags.reported).length;
  const included = hits.filter((h) => h.flags & HitFlags.included).length;
  const domZ = options.domZ ?? reported;
  if (!cutoffs) {
    for (const hit of hits) {
      const hitIncluded = (hit.flags & HitFlags.included) !== 0;
      for (const d of hit.domains) {
        d.reported = false;
        d.included = false;
        if (!(hit.flags & HitFlags.reported)) continue;
        d.reported = options.domScore === undefined ? d.p * domZ <= domE : d.bitscore >= options.domScore;
        d.included = hitIncluded &&
          (options.incdomScore === undefined ? d.p * domZ <= incdomE : d.bitscore >= options.incdomScore);
      }
    }
  }
  for (const hit of hits) {
    for (const d of hit.domains) d.cEvalue = d.p * domZ;
    // workaround_bug_h74(): keep only the best of identical domain alignments.
    if (hit.overlaps === 0) continue;
    for (let i = 0; i < hit.domains.length; i++) {
      for (let j = i + 1; j < hit.domains.length; j++) {
        const [a, b] = [hit.domains[i], hit.domains[j]];
        if (a.aliFrom !== b.aliFrom || a.aliTo !== b.aliTo) continue;
        const removed = a.bitscore >= b.bitscore ? b : a;
        removed.reported = false;
        removed.included = false;
      }
    }
  }
  // p7_tophits_SortBySortkey: sort key (-lnP when inclusion is by E-value, else score), then name.
  const byE = options.incScore === undefined && !cutoffs;
  const key = (h: Hit) => byE ? -h.lnP : h.score;
  const name = (h: Hit) => names ? names[h.target] : String(h.target);
  hits.sort((a, b) => key(b) - key(a) || compareBytes(name(a), name(b)));
  return {hits, z, domZ, reported, included};
}

/** C `strcmp` on UTF-8 bytes. */
function compareBytes(a: string, b: string): number {
  if (a === b) return 0;
  const [x, y] = [new TextEncoder().encode(a), new TextEncoder().encode(b)];
  for (let i = 0; i < Math.min(x.length, y.length); i++)
    if (x[i] !== y[i]) return x[i] - y[i];
  return x.length - y.length;
}
