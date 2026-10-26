import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {Fingerprint} from '../utils/chem-common';
import {SubstructureSearchType, getSearchQueryAndType} from '../constants';
import {CruxSegment, cruxReadsCxExtension, hashStrings} from './crux-service';
import {getCruxService, runCruxSearch} from './crux-searches';
import {SearchTypeKeys, SegmentSearchResult, searchTypeHits} from './jobs/search-types';

/** A search call's time per row of its segment (the first builds the segment's keys: two SMILES and a pattern
 * fingerprint per molecule, ~0.3 ms on Crux). */
const MS_PER_ROW = 2;

/** What the main thread keeps of a segment whose keys a worker holds: the rows Crux reads no molecule from (or whose
 * CXSMILES extension it does not read), which RDKit answers for, and their keys on RDKit (no RDKit objects); the
 * column version they were checked at and the hash of the values they were built from, extensions included (the
 * index's segments leave them out: an edit of an extension alone keeps the segment). */
interface SegmentState {
  unread: number[];
  rdkitKeys: SearchTypeKeys;
  version: number;
  hash: number;
}
const segmentStates = new WeakMap<CruxSegment, SegmentState>();

/** The column's values as the search types read them: whole (a CXSMILES with its extension), '' for an empty cell. */
function rawValues(col: DG.Column): string[] {
  return Array.from(col.toList() as (string | null)[], (v) => v ?? '');
}

/**
 * Crux counterpart of chemSubstructureSearchLibrary for Exact, Stereo agnostic, Included in and Not included in, with
 * the same contract (crux-searches.ts runCruxSearch): the column's index segments searched in their workers
 * (jobs/search-types.ts, each segment's keys built by its first search), RDKit's answer for the rows Crux cannot read
 * (the same code on the main thread's RDKit module), and RDKit's search for a search that fails. `includeMask` limits
 * Exact, Stereo agnostic and Included in, as RDKit's pattern screen does.
 */
export async function cruxSearchByType(col: DG.Column, molString: string, molBlockFailover: string,
  awaitAll: boolean, searchType: SubstructureSearchType, includeMask: BitArray | null, isSuperseded: () => boolean,
  rdkitSearch: (result: BitArray) => Promise<BitArray>): Promise<BitArray> {
  const searchKey = getSearchQueryAndType(molBlockFailover, searchType, Fingerprint.Morgan, 0);
  const masked = includeMask && searchType !== SubstructureSearchType.NOT_INCLUDED_IN ? includeMask : null;
  return runCruxSearch(col, searchKey, awaitAll, molString === '', isSuperseded, rdkitSearch, async (run) => {
    if (molString === '') // every row matches the empty query, as in RDKit's search
      return {search: Promise.resolve()};
    const service = await getCruxService();
    run.useService(service);
    const version = col.version;
    const index = await service.getIndex(col);
    let values: string[] | null = null; // read once a segment's keys are to be sent or checked
    let processed = 0;
    const search = service.eachSegment(index, async (segment, call) => {
      let sent: string[] | null = null;
      const send = (): string[] => {
        values ??= rawValues(col);
        // a CXSMILES whose extension Crux does not read as RDKit does: RDKit's answer (its row sent empty)
        return sent ??= values.slice(segment.start, segment.end).map((v) => cruxReadsCxExtension(v) ? v : '');
      };
      const job = (molecules: string[] | null) => call<SegmentSearchResult>('searchTypes::search',
        [segment.key, molecules, molString, molBlockFailover, searchType], MS_PER_ROW);
      let state = segmentStates.get(segment);
      if (state && state.version !== version) { // an edit: the keys of a segment whose values changed are built again
        values ??= rawValues(col);
        if (hashStrings(values, segment.start, segment.end) === state.hash)
          state.version = version;
        else
          state = undefined;
      }
      let answer = await job(state ? null : send());
      if (!answer) // its worker holds no keys of the segment (a worker replaced): sent again
        answer = await job(send());
      if (!answer)
        throw new Error('a Crux worker holds no search keys of a segment it was sent');
      if (!state) {
        const molecules = send();
        const unread = [...Array.from(answer.unread), ...molecules.flatMap((m, i) => !m && values![segment.start + i] ?
          [i] : [])].sort((a, b) => a - b);
        state = {unread, rdkitKeys: new SearchTypeKeys(unread.map((i) => values![segment.start + i])), version,
          hash: hashStrings(values!, segment.start, segment.end)};
        segmentStates.set(segment, state);
      }
      if (run.isCancelled())
        return;
      setHits(run.result, segment.start, answer.hits, masked);
      if (state.unread.length) {
        const rdkitHits = rdkitAnswer(state.rdkitKeys, molString, molBlockFailover, searchType);
        setHits(run.result, segment.start, Uint32Array.from(rdkitHits, (k) => state!.unread[k]), masked);
      }
      processed += segment.end - segment.start;
      run.progress(processed / col.length);
    }, run.isCancelled);
    return {search};
  });
}

/** RDKit's hits among the rows Crux cannot read; none where RDKit reads no query (Crux read it: the search goes on). */
function rdkitAnswer(keys: SearchTypeKeys, molString: string, molBlockFailover: string,
  searchType: SubstructureSearchType): number[] {
  try {
    return searchTypeHits(getRdKitModule(), keys, molString, molBlockFailover, searchType);
  } catch (e: any) {
    console.warn(`Chem | RDKit reads no ${searchType} query from ${molString}: no hits among the rows Crux cannot ` +
      `read (${e?.message ?? e})`);
    return [];
  }
}

function setHits(result: BitArray, start: number, hits: Uint32Array, includeMask: BitArray | null): void {
  for (const hit of hits) {
    if (!includeMask || includeMask.getBit(start + hit))
      result.setFast(start + hit, true);
  }
}

/**
 * Crux mode's Similar: the rows whose fingerprint of type `fp` (`fingerprints`: the column's, by the engine in effect,
 * cached per engine) has a Tanimoto similarity of at least `cutoff` to the query's (`queryFingerprint`), as RDKit's
 * search computes it (the same Tanimoto, no index); with the contract of chemSubstructureSearchLibrary
 * (crux-searches.ts runCruxSearch), RDKit's search for a search that fails.
 */
export async function cruxSimilaritySearch(col: DG.Column, molString: string, molBlockFailover: string,
  awaitAll: boolean, cutoff: number, fp: Fingerprint, isSuperseded: () => boolean,
  rdkitSearch: (result: BitArray) => Promise<BitArray>, fingerprints: () => Promise<(BitArray | null)[]>,
  queryFingerprint: () => BitArray): Promise<BitArray> {
  const searchKey = getSearchQueryAndType(molBlockFailover, SubstructureSearchType.IS_SIMILAR, fp, cutoff);
  return runCruxSearch(col, searchKey, awaitAll, molString === '', isSuperseded, rdkitSearch, async (run) => {
    if (molString === '') // every row matches the empty query, as in RDKit's search
      return {search: Promise.resolve()};
    const query = queryFingerprint(); // throws for a query no engine reads, as RDKit's search does
    const search = fingerprints().then((fps) => {
      if (run.isCancelled())
        return;
      for (let i = 0; i < fps.length; i++) {
        if (fps[i] && tanimotoSimilarity(fps[i]!, query) >= cutoff)
          run.result.setFast(i, true);
      }
    });
    return {search};
  });
}
