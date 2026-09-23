import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {CruxDataset, CruxPool} from '../../wasm/crux/pool.js';
import {getSearchProgressEventName, getSearchQueryAndType, getTerminateEventName,
  SubstructureSearchType} from '../constants';
import {Fingerprint} from '../utils/chem-common';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {PackageFunctions} from '../package';
import {chemSubstructureSearchLibrary, FILTER_TYPES, subscribeToColumnChanges} from '../chem-searches';

const PROGRESS_MILESTONES = [1, 5, 10, 20, 40, 70]; // Completion is the seventh and final update.
const MIN_PROGRESS_INTERVAL_MS = 150;

/** A filter owns its pool; indexes survive sketch edits and workers stop when the filter detaches. */
export class CruxSubstructureService {
  private pool?: Promise<CruxPool>;
  private cache?: {column: DG.Column; version: number; length: number; dataset: Promise<CruxDataset>};
  private buildTail: Promise<void> = Promise.resolve();
  private lifetime = new AbortController();
  private query?: AbortController;

  private getPool(): Promise<CruxPool> {
    this.lifetime.signal.throwIfAborted();
    const workers = Math.max(1, Math.min(8, Math.floor((navigator.hardwareConcurrency || 2) / 2)));
    this.pool ??= CruxPool.create({workers})
      .then(async (pool) => {
        if (this.lifetime.signal.aborted) {
          await pool.dispose();
          this.lifetime.signal.throwIfAborted();
        }
        return pool;
      }).catch((error) => {
        this.pool = undefined;
        throw error;
      });
    return this.pool;
  }

  private getDataset(column: DG.Column): Promise<CruxDataset> {
    const version = column.version;
    const length = column.length;
    if (this.cache && this.cache.column.dart === column.dart &&
      this.cache.version === version && this.cache.length === length)
      return this.cache.dataset;

    // Snapshot before yielding: row positions must correspond to the version used as the cache key.
    const categories = column.categories.slice();
    const rowCategories = column.getRawData().slice(0, length);
    const dataset = this.buildTail.then(async () => {
      const pool = await this.getPool();
      let values = categories;
      if (categories.some((mol) => DG.chem.isMolBlock(mol)))
        values = await (await getRdKitService()).convertMolNotation(categories, DG.chem.Notation.Smiles);
      this.lifetime.signal.throwIfAborted();
      // The worker transport is newline delimited; malformed cells must not shift subsequent row indexes.
      values = values.map((mol) => mol && !/[\r\n]/.test(mol) ? mol : '');
      const smiles = new Array<string>(length);
      for (let i = 0; i < length; i++)
        smiles[i] = values[rowCategories[i]] ?? '';
      return pool.load(smiles, {buildIndex: true, shardSize: 25_000, chunkSize: 25_000,
        signal: this.lifetime.signal});
    });
    this.cache = {column, version, length, dataset};
    this.buildTail = dataset.then(() => {}, () => {
      if (this.cache?.dataset === dataset)
        this.cache = undefined;
    });
    return dataset;
  }

  /** Uses the same mutable bit array and progress/completion events as the RDKit filter endpoint. */
  async search(column: DG.Column, molecule: string, smarts: string, fp: Fingerprint, cutoff: number,
    awaitAll = false, includeMask?: BitArray): Promise<BitArray> {
    this.query?.abort();
    this.lifetime.signal.throwIfAborted();
    subscribeToColumnChanges(column);
    const queryMol = getQueryMolSafe(molecule, smarts, PackageFunctions.getRdKitModule());
    if (!queryMol) {
      grok.shell.error('Crux: Search pattern cannot be set');
      throw new Error('Crux: Search pattern cannot be set');
    }
    let querySmarts: string;
    let hasRadicals: boolean;
    try {
      querySmarts = queryMol.get_smarts();
      // RDKit's mol-to-SMARTS serialization omits radical electron constraints.
      const json = JSON.parse(queryMol.get_json());
      hasRadicals = json.molecules.some((mol: {atoms: {nRad?: number}[]}) =>
        mol.atoms.some((atom) => (atom.nRad ?? 0) > 0));
    } finally {
      queryMol.delete();
    }
    if (hasRadicals) {
      return chemSubstructureSearchLibrary(column, molecule, smarts, FILTER_TYPES.substructure,
        false, awaitAll, SubstructureSearchType.CONTAINS, cutoff, fp, includeMask);
    }
    const controller = new AbortController();
    this.query = controller;
    const signal = controller.signal;
    const matches = new BitArray(column.length);
    const tableName = column.dataFrame?.name ?? '';
    const progressEvent = getSearchProgressEventName(tableName, column.name);
    const terminateEvent = getTerminateEventName(tableName, column.name);
    const queryId = getSearchQueryAndType(smarts, SubstructureSearchType.CONTAINS, fp, cutoff);
    const termination = grok.events.onCustomEvent(terminateEvent).subscribe((id: string) => {
      if (id === queryId)
        controller.abort();
    });

    const run = async () => {
      let milestone = 0;
      let lastUpdateTime = -Infinity;
      try {
        // Give the filter time to subscribe, including when the index is already cached.
        await new Promise<void>((resolve) => setTimeout(resolve, 0));
        signal.throwIfAborted();
        const dataset = await this.getDataset(column);
        signal.throwIfAborted();
        for await (const batch of dataset.substructureSearchStream(querySmarts, {signal})) {
          signal.throwIfAborted();
          for (const index of batch.indices) {
            if (!includeMask || includeMask.getBit(index))
              matches.setBit(index, true);
          }
          const progress = 100 * batch.shardsDone / batch.shardCount;
          const now = performance.now();
          // A progress event also rebuilds the bitset and refilters/redraws the table.
          if (!awaitAll && progress < 100 && milestone < PROGRESS_MILESTONES.length &&
            progress >= PROGRESS_MILESTONES[milestone] && now - lastUpdateTime >= MIN_PROGRESS_INTERVAL_MS) {
            while (milestone < PROGRESS_MILESTONES.length && progress >= PROGRESS_MILESTONES[milestone])
              milestone++;
            lastUpdateTime = now;
            grok.events.fireCustomEvent(progressEvent, progress);
          }
        }
        signal.throwIfAborted();
      } catch (error) {
        if (!signal.aborted && !this.lifetime.signal.aborted)
          grok.shell.error(`Crux substructure search failed: ${error instanceof Error ? error.message : error}`);
        if (awaitAll)
          throw error;
      } finally {
        termination.unsubscribe();
        if (!awaitAll && !signal.aborted && !this.lifetime.signal.aborted) {
          grok.events.fireCustomEvent(progressEvent, 100);
          grok.events.fireCustomEvent(terminateEvent, queryId);
        }
      }
    };
    if (awaitAll)
      await run();
    else
      void run();
    return matches;
  }

  dispose(): void {
    this.query?.abort();
    this.lifetime.abort();
    void this.pool?.then((pool) => pool.dispose()).catch(() => {});
    this.cache = undefined;
  }
}
