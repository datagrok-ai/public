import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_package} from '../package';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {Fingerprint} from '../utils/chem-common';
import {SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType,
  getTerminateEventName} from '../constants';
import {CruxColumnIndex, CruxSegment, CruxService} from './crux-service';
import {isCruxRequested} from '../engine/chem-engine';

export enum SubstructureSearchEngine {
  RDKit = 'RDKit',
  Crux = 'Crux',
}

/** At most this many filter updates before the final one. */
const MAX_INTERMEDIATE_UPDATES = 6;
/** Intermediate updates keep this far from the search start, from each other and from the expected end. */
const MIN_UPDATE_GAP_MS = 100;

let engineOverride: SubstructureSearchEngine | null = null;
let cruxService: Promise<CruxService> | null = null;

/** Tests: overrides the Chemical Engine's choice of the search engine until reset with null. */
export function setSubstructureSearchEngine(engine: SubstructureSearchEngine | null): void {
  engineOverride = engine;
}

/**
 * The engine of Contains / Not contains: the Chemical Engine's (Crux in Crux mode, RDKit in RDKit mode, or once the
 * Crux module has failed to load); it has no setting of its own since 2026-10-22, and a value a server saved for the
 * retired SubstructureSearchEngine is not read.
 */
export function getSubstructureSearchEngine(): SubstructureSearchEngine {
  if (engineOverride)
    return engineOverride;
  return isCruxRequested() ? SubstructureSearchEngine.Crux : SubstructureSearchEngine.RDKit;
}

export function getCruxService(): Promise<CruxService> {
  cruxService ??= CruxService.create().catch((e) => {
    cruxService = null;
    throw e;
  });
  return cruxService;
}

/** Tests: the next search starts the Crux service anew (from the loader); the current one's workers are let go. */
export async function resetCruxServiceForTests(): Promise<void> {
  const service = await cruxService?.catch(() => null);
  cruxService = null;
  service?.terminateForTests();
}

export function isCruxSearch(searchType: SubstructureSearchType): boolean {
  return (searchType === SubstructureSearchType.CONTAINS || searchType === SubstructureSearchType.NOT_CONTAINS) &&
    getSubstructureSearchEngine() === SubstructureSearchEngine.Crux;
}

/**
 * Whether Crux searches for the query: the empty one (the pre-search on filter attach, which builds the index), or one
 * Crux reads (getQueryMolSafe gives a molecule on Chem's Crux module), which Crux then matches as RDKit does. Otherwise
 * (a query Crux does not read, Crux not available) RDKit searches.
 */
export async function cruxReadsQuery(molString: string, molBlockFailover: string): Promise<boolean> {
  try {
    const service = await getCruxService();
    if (molString === '' || service.readsQuery(molString, molBlockFailover))
      return true;
    _package.logger.debug(`Crux reads no query from ${molString}, searching with RDKit`);
    return false;
  } catch (e: any) {
    console.warn(`Chem | crux search is not available, searching with RDKit: ${e?.message ?? e}`);
    return false;
  }
}

/**
 * Crux counterpart of chemSubstructureSearchLibrary (Contains / Not contains) with the same contract: the returned
 * bit array fills in as segments are searched, reported by the progress event (at most 7 times per search, the last
 * one at 100%), and the terminate event, keyed by query, both ends the search and cancels it. A search that is not
 * awaited stops, announcing its end, once `isSuperseded` tells a newer search on the column started. If crux fails
 * (a worker out of memory, say), crux restarts, which frees its memory, and `rdkitSearch` fills the bit array instead.
 */
export async function cruxSubstructureSearch(col: DG.Column, molString: string,
  molBlockFailover: string, awaitAll: boolean, searchType: SubstructureSearchType,
  includeMask: BitArray | null, isSuperseded: () => boolean,
  rdkitSearch: (result: BitArray) => Promise<BitArray>): Promise<BitArray> {
  const dfName = col.dataFrame?.name ?? '';
  const terminateEventName = getTerminateEventName(dfName, col.name);
  const progressEventName = getSearchProgressEventName(dfName, col.name);
  const searchKey = getSearchQueryAndType(molBlockFailover, searchType, Fingerprint.Morgan, 0);
  const result = new BitArray(col.length);
  // the pre-search on filter attach (an empty query): every row matches, the index gets built in the background
  const preSearch = molString === '';
  if (preSearch)
    result.setAll(true);
  let terminated = false;
  let failed = false;
  const terminateSub = awaitAll ? null : grok.events.onCustomEvent(terminateEventName).subscribe((key: string) => {
    terminated ||= key === searchKey;
  });
  const superseded = () => !awaitAll && !preSearch && isSuperseded();
  const isCancelled = () => terminated || failed || superseded();

  // molecules crux cannot parse get the verdict of the RDKit search, which parses more leniently
  let rdkitQuery: RDMol | null | undefined;
  const rdkitMatch = (molecule: string): boolean | null => {
    const mol = molecule ? getMolSafe(molecule, {}, getRdKitModule()).mol : null;
    try {
      rdkitQuery ??= mol ? getQueryMolSafe(molString, molBlockFailover, getRdKitModule()) : undefined;
      return mol ? !!rdkitQuery && mol.get_substruct_match(rdkitQuery) !== '{}' : null;
    } finally {
      mol?.delete();
    }
  };

  let service: CruxService | null = null;
  let generation = 0;
  /** Runs the search on RDKit; `restart` when the crux workers failed, which frees their memory. */
  const fallBack = async (e: any, restart: boolean): Promise<BitArray> => {
    rdkitQuery?.delete();
    failed = true;
    if (restart)
      service?.restart(generation);
    if (terminated || superseded() || preSearch) {
      terminateSub?.unsubscribe();
      if (!terminated)
        grok.events.fireCustomEvent(terminateEventName, searchKey);
      return result;
    }
    console.warn(`Chem | crux substructure search failed, searching with RDKit: ${e?.message ?? e}`);
    result.setAll(false);
    try {
      return await rdkitSearch(result);
    } finally {
      // a terminate that came while the RDKit search was starting did not reach it
      terminateSub?.unsubscribe();
      if (terminated)
        grok.events.fireCustomEvent(terminateEventName, searchKey);
    }
  };

  let index: CruxColumnIndex;
  let search: Promise<void>;
  try {
    service = await getCruxService();
    generation = service.generation;
    index = await service.getIndex(col);
    if (preSearch)
      search = service.build(index, isCancelled);
    else {
      const updates = new FilterUpdateScheduler((fraction) =>
        grok.events.fireCustomEvent(progressEventName, fraction * 100));
      let processed = 0;
      search = service.search(index, molString, molBlockFailover, (segment, hits) => {
        setMatches(result, segment, hits, searchType, includeMask, rdkitMatch);
        processed += segment.end - segment.start;
        if (!awaitAll)
          updates.onProgress(processed / col.length);
      }, isCancelled);
    }
  } catch (e: any) {
    return fallBack(e, false);
  }

  const finished = search.then(() => {
    terminateSub?.unsubscribe();
    rdkitQuery?.delete();
    if (terminated)
      return;
    if (!superseded()) {
      reportUnparsed(index, col);
      if (!awaitAll && !preSearch)
        grok.events.fireCustomEvent(progressEventName, 100);
    }
    grok.events.fireCustomEvent(terminateEventName, searchKey);
  }, (e) => fallBack(e, true));
  if (awaitAll)
    await finished;
  else
    finished.catch(() => {});
  return result;
}

function setMatches(result: BitArray, segment: CruxSegment, hits: Uint32Array, searchType: SubstructureSearchType,
  includeMask: BitArray | null, rdkitMatch: (molecule: string) => boolean | null): void {
  // per segment row: 1 - matches, 2 - does not match, 0 - not a molecule
  const verdicts = new Uint8Array(segment.end - segment.start).fill(2);
  for (const hit of hits)
    verdicts[hit] = 1;
  for (let k = 0; k < segment.failed!.length; k++) {
    const match = rdkitMatch(segment.unparsed![k]);
    if (match === null)
      segment.unparsed![k] = '';
    verdicts[segment.failed![k]] = match === null ? 0 : match ? 1 : 2;
  }
  const wanted = searchType === SubstructureSearchType.CONTAINS ? 1 : 2;
  for (let i = 0; i < verdicts.length; i++) {
    const row = segment.start + i;
    if (verdicts[i] === wanted && (!includeMask || includeMask.getBit(row)))
      result.setFast(row, true);
  }
}

const reportedIndexes = new WeakSet<CruxColumnIndex>();

/** Tells once per column index which molecules crux could not parse (RDKit decides for them). */
function reportUnparsed(index: CruxColumnIndex, col: DG.Column): void {
  if (reportedIndexes.has(index) || index.segments.some((s) => !s.unparsed))
    return;
  reportedIndexes.add(index);
  const unparsed = index.segments.flatMap((s) => s.unparsed!).filter((m) => m !== '');
  if (unparsed.length > 0) {
    console.info(`Chem | crux could not parse ${unparsed.length} molecules of ${col.name}, RDKit decides for them: ` +
      unparsed.slice(0, 5).join(' '));
  }
}

/**
 * Spreads the intermediate filter updates of a search: the first one as soon as results come in (a quick search
 * shows them once, at the end), the next ones at geometrically growing progress.
 */
class FilterUpdateScheduler {
  private sent = 0;
  private next = 0;
  private ratio = 2;
  private readonly start = performance.now();
  private last = this.start;

  constructor(private readonly update: (fraction: number) => void) {}

  onProgress(fraction: number): void {
    if (this.sent >= MAX_INTERMEDIATE_UPDATES || fraction <= 0 || fraction >= 1 || fraction < this.next)
      return;
    const now = performance.now();
    if (now - this.last < MIN_UPDATE_GAP_MS || (now - this.start) * (1 - fraction) / fraction < MIN_UPDATE_GAP_MS)
      return;
    if (this.sent === 0)
      this.ratio = Math.max(2, Math.pow(1 / fraction, 1 / MAX_INTERMEDIATE_UPDATES));
    this.update(fraction);
    this.sent++;
    this.last = now;
    this.next = fraction * this.ratio;
  }
}
