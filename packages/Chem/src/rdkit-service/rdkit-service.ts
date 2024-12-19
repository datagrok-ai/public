import * as DG from 'datagrok-api/dg';
import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import {Fingerprint, rdKitFingerprintToBitArray} from '../utils/chem-common';
import {RuleId} from '../panels/structural-alerts';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {IFpResult} from './rdkit-service-worker-similarity';
import {LockedEntity} from '../utils/locked-entitie';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {getRdKitModule} from '../package';
import {SubstructureSearchType} from '../constants';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {getRDKitFpAsUint8Array} from '../chem-searches';
import {IMmpFragmentsResult, InverseSubstructureRes, IRGroupAnalysisResult} from './rdkit-service-worker-substructure';
import { ISubstruct } from '@datagrok-libraries/chem-meta/src/types';
export interface IParallelBatchesRes {
  getProgress: () => number,
  setTerminateFlag: () => void,
  getTerminateFlag: () => boolean,
  promises: Promise<void>[]
}

export type SubstructureSearchWithFpResult = {
  bitArray: BitArray,
  fpsRes: IFpResult | null,
};

export type SubstructureSearchBatchResult = {
  matches: BitArray;
  fpRes: IFpResult | null;
  fpCreated: boolean
};

export class RdKitService {
  workerCount: number;
  _initWaiters?: Promise<any>[];
  timesInitialized = 0;
  parallelWorkers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;
  moleculesSegmentsLengths: Uint32Array;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this.workerCount = Math.max(1, cpuLogicalCores - 2);
    this.moleculesSegmentsLengths = new Uint32Array(this.workerCount);
  }

  async init(webRoot: string): Promise<void> {
    if (!this._initWaiters) {
      this._initWaiters = [];
      for (let i = 0; i < this.workerCount; ++i) {
        const workerClient = new RdKitServiceWorkerClient();
        this.parallelWorkers[i] = workerClient;
        this._initWaiters.push(workerClient.moduleInit(webRoot));
      }
    }
    await Promise.all(this._initWaiters);
    if (this.timesInitialized++ === 0)
      console.log('RDKit Service was initialized');
  }

  /**
   * Performs parallelization of function execution using web-workers.
   * @async
   * @param {function (workerIdx: number, workerCount: number): Promise<TMap>} map - splits the data
   * by number of workers and calls function from rdkit service worker client (basicaly action which we need
   * to perform inside worker - getFingerprints, searchSubstructure etc.)
   * @param {function (_: TMap[]): TReduce} reduce - fucntion which combines results collected from web workers
   * into single result
   * */
  async _doParallel<TMap, TReduce>(
    map: (workerIdx: number, workerCount: number) => Promise<TMap>,
    reduce: (_: TMap[]) => TReduce = (_: TMap[]) => [] as TReduce): Promise<TReduce> {
    const promises = [];
    const workerCount = this.workerCount;
    for (let workerIdx = 0; workerIdx < workerCount; workerIdx++)
      promises[workerIdx] = map(workerIdx, workerCount);

    const data = await Promise.all(promises);
    return reduce(data);
  }


  async _doParallelBatches<TData, TRes, BatchRes>(
    data: TData[],
    res: TRes,
    updateRes: (batchRes: BatchRes, res: TRes, length: number, index: number) => void,
    map: (batch: TData[], workerIdx: number, workerCount: number, batchStartIdx: number) => Promise<BatchRes>,
    progressFunc: (progress: number) => void): Promise<IParallelBatchesRes> {
    let terminateFlag = false;

    const setTerminateFlag = () => {terminateFlag = true;};
    const getTerminateFlag = () => {return terminateFlag;};
    const t = this;
    const dataLength = data.length;
    const increment = Math.floor(Math.max(500 / this.workerCount, 20));
    const incrementMultiplier = 1.05; // each iteration increment is multiplied by this value to increase sub_batch size
    const workingIndexes = new Array<{start: number, increment: number}>(this.workerCount)
      .fill({start: 0, increment});
    // const sizePerWorker = Math.floor(dataLength / this.workerCount);

    /*distribute data between workers, makes sure that same molecules always end up in the same worker,
    important for caching */
    for (let i = 0; i < this.workerCount; i++) {
      workingIndexes[i] = {
        start: i * increment,
        increment,
      }; // if we have 3 workers per se, the working indexes are disted in this way between workers:
      // | 1 | 2 | 3 |  1  |  2  |  3  |   1   |   2   |   3   |    1    |    2    |    3    |
    }

    const lockedCounter = new LockedEntity(0);
    // array that stores how many iterations each worker did
    const workerIterationCounter = new Array<number>(this.workerCount).fill(0);
    let iterationCheckPoint = 1; // next checkpoint of iterations for sending progress
    let iterationsPerProgress = 1; // how many iterations should be done before sending progress
    let processedMolecules = 0;

    // array with worker indexes containing map of iteration number to updateRes closure
    const updateResMapArray = new Array(this.workerCount).fill(0)
      .map(() => new Map<number, {batchResult: BatchRes, curWorkerStart: number, l: number}>());
    const promises = t.parallelWorkers.map((_, idx) => {
      const post = async () => {
        if (workingIndexes[idx].start >= dataLength || terminateFlag)
          return;

        const part = data.slice(workingIndexes[idx].start,
          Math.min(dataLength, workingIndexes[idx].start + workingIndexes[idx].increment));
        const batchResult = await map(part, idx, t.parallelWorkers.length, workingIndexes[idx].start);
        if (terminateFlag)
          return;

        await lockedCounter.unlockPromise();
        lockedCounter.lock();
        const curWorkerStart = workingIndexes[idx].start;
        updateResMapArray[idx].set(workerIterationCounter[idx],
          {batchResult: batchResult, curWorkerStart: curWorkerStart, l: part.length});

        workingIndexes[idx].start += // worker with index idx will have this.workerCount - idx chunks in front with size
          part.length * (this.workerCount - idx) + // of old increment and idx amount of chunks with new increment size
          idx * (Math.floor(part.length * (incrementMultiplier))); // hence the two terms
        processedMolecules = lockedCounter.value;
        const end = Math.min(processedMolecules + workingIndexes[idx].increment, dataLength);
        workerIterationCounter[idx] += 1;
        const updatePromises: Promise<void>[] = [];
        if (workerIterationCounter.every((it) => it >= iterationCheckPoint)) {
          // update all results:
          for (let i = 0; i < this.workerCount; i++) {
            const updateResMap = updateResMapArray[i];
            for (const it of updateResMap.keys()) {
              if (it <= iterationCheckPoint) {
                const br = updateResMap.get(it)!;
                updatePromises.push(new Promise<void>((resolveUpdate) => {
                  setTimeout(() => {
                    updateRes(br.batchResult, res, br.l, br.curWorkerStart);
                    updateResMap.delete(it);
                    resolveUpdate();
                  });
                }));
              }
            }
          }
          iterationsPerProgress *= 1.05; // increase iterations per progress
          iterationCheckPoint += iterationsPerProgress;
          iterationCheckPoint = Math.floor(iterationCheckPoint);
          const progressFraction = processedMolecules/dataLength;
          Promise.all(updatePromises).then(() => progressFunc(progressFraction));
        }

        workingIndexes[idx].increment = Math.floor(workingIndexes[idx].increment * incrementMultiplier);
        lockedCounter.value = end;
        lockedCounter.release();
        // after the locked counter is released, wait for all update promises to finish;
        await Promise.all(updatePromises);
        await post();
      };
      return post();
    });
    // some workers might have finished before others, so we need to update the results with what is left
    Promise.all(promises).then(() => {
      updateResMapArray.forEach((leftMap) => {
        for (const it of leftMap.keys()) {
          const br = leftMap.get(it)!;
          updateRes(br.batchResult, res, br.l, br.curWorkerStart);
          leftMap.delete(it);
        }
      });
      progressFunc(1);
    });
    const getProgress = () => processedMolecules / dataLength;
    return {getProgress, setTerminateFlag, getTerminateFlag, promises};
  }


  /**
   * Calls _doParallel with pre-defined map function which splits data by number of workers
   * @async
   * @param {Array<T>} molecules - list of molecules to split by workers
   * @param {function (workerIdx: number, workerCount: number): Promise<TMap>} workerFunc - function
   * from rdkit service worker client (basicaly action which we need to perform inside worker - getFingerprints,
   * searchSubstructure etc.)
   * @param {function (_: TMap[]): TReduce} reduce - function which combines results collected from web workers
   * into single result
   * */
  async _initParallelWorkers<TMap, TReduce, T>(molecules: Array<T>,
    workerFunc: (workerIdx: number, moleculesSegment: Array<T>) => Promise<TMap>,
    reduce: (_: TMap[]) => TReduce): Promise<TReduce> {
    const t = this;
    return this._doParallel(
      (workerIdx: number, nWorkers: number) => {
        const length = molecules.length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segment = workerIdx < (nWorkers - 1) ?
          molecules.slice(workerIdx * segmentLength, (workerIdx + 1) * segmentLength) :
          molecules.slice(workerIdx * segmentLength, length);
        t.moleculesSegmentsLengths![workerIdx] = segment.length;
        return workerFunc(workerIdx, segment);
      },
      reduce,
    );
  }

  /**
   * Calls _doParallel with pre-defined map function which splits data by number of workers
   * @async
   * @param {Array<Array<T>>} data - list of molecules to split by workers
   * @param {function (workerIdx: number, workerCount: number): Promise<TMap>} workerFunc - function
   * from rdkit service worker client (basicaly action which we need to perform inside worker - getFingerprints,
   * searchSubstructure etc.)
   * @param {function (_: TMap[]): TReduce} reduce - function which combines results collected from web workers
   * into single result
   * */
  async _initParallelWorkersArray<TMap, TReduce, T>(data: Array<Array<T>>,
    workerFunc: (workerIdx: number, dataSegment: Array<Array<T>>) => Promise<TMap>,
    reduce: (_: TMap[]) => TReduce): Promise<TReduce> {
    const t = this;
    const lengthAll = data.length;
    return this._doParallel(
      (workerIdx: number, nWorkers: number) => {
        const length = data[0].length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segmentArray = Array<Array<T>>(lengthAll);
        for (let i = 0; i < lengthAll; i++) {
          const segment = workerIdx < (nWorkers - 1) ?
            data[i].slice(workerIdx * segmentLength, (workerIdx + 1) * segmentLength) :
            data[i].slice(workerIdx * segmentLength, length);
          segmentArray[i] = segment;
        }

        t.moleculesSegmentsLengths![workerIdx] = segmentArray[0].length;
        return workerFunc(workerIdx, segmentArray);
      },
      reduce,
    );
  }

  /**
   * Fills array of mols in each worker
   * @async
   * @param {string[]} molecules - list of molecules to save in each worker
   * */
  async initMoleculesStructures<TReduce>(molecules: string[])
    : Promise<TReduce | void> {
    return this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
      this.parallelWorkers[i].initMoleculesStructures(segment),
    () => {});
  }

  /**
   * Filters molecules by substructure
   * @async
   * @param {string} query - smiles/molblock to filter by
   * @param {string} queryMolBlockFailover - smart to filter by (is used if creation of RDMol object from query
   * parameter failed)
   * @param {SubstructureSearchWithFpResult} result - object where to put the results of substructure search and
   * created fp
   * @param {function (progress: number): void} progressFunc - function which should be run on each progress iteration
   * @param {boolean} molecules - list of molecules to search
   * @param {boolean} createSmiles - whether canonical smiles should be created besides fp (should be true if molecules
   * are not canonical smiles)
   * */
  async searchSubstructureWithFps(query: string, queryMolBlockFailover: string, result: SubstructureSearchWithFpResult,
    progressFunc: (progress: number) => void, molecules: string[], createSmiles = false,
    searchType = SubstructureSearchType.CONTAINS, simCutOff = 0.8, fp = Fingerprint.Morgan, afterBatchCalculated = () => {}) {
    const queryMol = searchType === SubstructureSearchType.IS_SIMILAR ? getMolSafe(query, {}, getRdKitModule()).mol :
      getQueryMolSafe(query, queryMolBlockFailover, getRdKitModule());
    const fpType = searchType === SubstructureSearchType.IS_SIMILAR ? fp : Fingerprint.Pattern;
    if (!queryMol)
      throw new Error(`Chem | Invalid search pattern: ${query}`);
    let fpRdKit: Uint8Array;
    try {
      fpRdKit = searchType === SubstructureSearchType.IS_SIMILAR ? getRDKitFpAsUint8Array(queryMol, fpType):
        queryMol.get_pattern_fp_as_uint8array();
    } catch (e: any) {
      throw new Error(`Chem | Substructure Search failed with error: ${e.toString()}`);
    } finally {
      queryMol?.delete();
    }

    const updateRes = (batchRes: SubstructureSearchBatchResult, res: SubstructureSearchWithFpResult, length: number,
      index: number) => {
      if (!res.fpsRes) {
        res.fpsRes = {
          fps: new Array<Uint8Array | null>(molecules.length).fill(null),
          smiles: createSmiles ? new Array<string | null>(molecules.length).fill(null) : null,
        };
      }

      for (let j = 0; j < length; j++) {
        //bit = !!(batchRes!.matches[Math.floor(j / 32)] >> j % 32 & 1);
        res.bitArray.setBit(index + j, batchRes.matches.getBit(j));
        res.fpsRes.fps[index + j] = batchRes.fpRes!.fps[j];
        if (createSmiles) {
              res.fpsRes.smiles![index + j] = batchRes.fpRes!.smiles![j];
              batchRes.fpRes!.smiles![j] = null;
        }
            batchRes.fpRes!.fps[j] = null;
      }
      batchRes.fpRes = null;
      if (batchRes.fpCreated)
        afterBatchCalculated();
    };

    return this._doParallelBatches(molecules, result, updateRes, async (batch, workerIdx,
      _workerCount, batchStartIdx) => {
      let fpResult: IFpResult;
      let fpCreated = false;
      if (!result.fpsRes || !result.fpsRes.fps[batchStartIdx + batch.length - 1]) {
        fpResult = await this.parallelWorkers[workerIdx].getFingerprints(fpType, batch, createSmiles);
        fpCreated = true;
      }
      else {
        fpResult = {
          fps: result.fpsRes.fps.slice(batchStartIdx, batchStartIdx + batch.length),
          smiles: createSmiles ? result.fpsRes.smiles!.slice(batchStartIdx, batchStartIdx + batch.length) : batch,
        };
      }

      if (query === '')
        return {matches: new BitArray(batch.length, true), fpRes: fpResult, fpCreated: fpCreated};
      let finalBitArray = new BitArray(batch.length, false);
      if (searchType !== SubstructureSearchType.IS_SIMILAR) {
        let filteredMolecules: string[] = [];
        let patternFpFilterBitArray: BitArray | null = null;
        if (searchType !== SubstructureSearchType.NOT_CONTAINS && searchType !== SubstructureSearchType.NOT_INCLUDED_IN) {
          // *********** FILTERING using fingerprints
          patternFpFilterBitArray = this.filterByPatternFps(searchType, batch, fpRdKit, fpResult);
          filteredMolecules = this.filterMoleculesByBitArray(patternFpFilterBitArray, batch, fpResult, createSmiles);
        } else
          filteredMolecules = createSmiles ? fpResult.smiles! as string[] : batch;

        // *********** DONE FILTERING using fingerprints if necessary
        // filter using substruct search on already prefiltered dataset
        const substructRes: Uint32Array = await this.parallelWorkers[workerIdx]
          .searchSubstructure(query, queryMolBlockFailover, filteredMolecules!, searchType);

        const matchesBitArray = BitArray.fromUint32Array(filteredMolecules.length, substructRes);
        if (searchType !== SubstructureSearchType.NOT_CONTAINS && searchType !== SubstructureSearchType.NOT_INCLUDED_IN)
          this.restorePrefilteredMoleculesIndexes(patternFpFilterBitArray!, matchesBitArray, finalBitArray);
        else
          finalBitArray = matchesBitArray;
      } else {
        // ************* PERFORM SIMILARITY SEARCH
        try {
          checkEl:
          for (let i = 0; i < batch.length; ++i) {
            if (fpResult.fps[i]) {
              const simScore = tanimotoSimilarity(rdKitFingerprintToBitArray(fpResult.fps[i]!)!, rdKitFingerprintToBitArray(fpRdKit!)!);
              if (simScore >= simCutOff)
                finalBitArray.setBit(i, true);
            }
          }
        } catch (e: any) {
          throw new Error(`Chem | Similarity Search failed with error: ${e.toString()}`);
        }
      }

      return {
        matches: finalBitArray,
        fpRes: fpResult,
        fpCreated: fpCreated
      };
    }, progressFunc);
  }

  // restore the indexes of prefiltered molecules on the whole dataset
  restorePrefilteredMoleculesIndexes(patternFpFilterBitArray: BitArray, matchesBitArray: BitArray, finalBitArray: BitArray) {
    // restore the indexes of prefiltered molecules on the whole dataset
    let matchesCounter = 0;
    for (let i = -1; (i = patternFpFilterBitArray!.findNext(i)) != -1;) {
      if (matchesBitArray.getBit(matchesCounter))
        finalBitArray.setBit(i, true);
      matchesCounter++;
    }
  }

  filterByPatternFps(searchType: SubstructureSearchType, batch: string[], fpRdKit: Uint8Array,
    fpResult: IFpResult): BitArray {
    const superStructSearch = !!(searchType === SubstructureSearchType.INCLUDED_IN);
    const patternFpUint8Length = 256;
    const patternFpFilterBitArray = new BitArray(batch.length, false);
    checkEl:
    for (let i = 0; i < batch.length; ++i) {
      if (fpResult.fps[i]) {
        for (let j = 0; j < patternFpUint8Length; ++j) {
          const bitToCompare = superStructSearch ? fpResult.fps[i]![j] : fpRdKit[j];
          if ((fpResult.fps[i]![j] & fpRdKit[j]) != bitToCompare)
            continue checkEl;
        }
        patternFpFilterBitArray.setFast(i, true);
      }
    }
    return patternFpFilterBitArray;
  }

  filterMoleculesByBitArray(patternFpFilterBitArray: BitArray, batch: string[], fpResult: IFpResult, createSmiles?: boolean): string[] {
    const filteredMolecules = Array<string>(patternFpFilterBitArray.trueCount());
    let counter = 0;
    for (let i = -1; (i = patternFpFilterBitArray.findNext(i)) !== -1;) {
      filteredMolecules[counter] = createSmiles ? fpResult.smiles![i]! : batch[i];
      counter++;
    }
    return filteredMolecules;
  }


  /**
   * Returns fingerprints for provied molecules
   * @async
   * @param {Fingerprint} fingerprintType - type of fingerprint (Morgan or Pattern)
   * @param {boolean} molecules - optional parameter, if passed RDMols objects are created on the fly (array of
   * predefined RDMols is not used)
   * @param {boolean} getCanonicalSmiles - optional parameter, if passed, in addition to fps function returns
   * canonical smiles
   * */
  async getFingerprints(fingerprintType: Fingerprint, molecules?: string[],
    getCanonicalSmiles?: boolean): Promise<IFpResult> {
    const t = this;
    const getResult = (data: IFpResult[]): IFpResult => {
      return {
        fps: ([] as Array<Uint8Array | null>).concat(...data.map(((it) => it.fps))),
        smiles: getCanonicalSmiles ? ([] as Array<string | null>).concat(...data.map(((it) => it.smiles))) : null,
      };
    };
    const res = molecules ?
      await this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
        t.parallelWorkers[i].getFingerprints(fingerprintType, segment, getCanonicalSmiles),
      (data: IFpResult[]) => {
        return getResult(data);
      }) :
      await this._doParallel(
        (i: number, _: number) => {
          return t.parallelWorkers[i].getFingerprints(fingerprintType, molecules, getCanonicalSmiles);
        },
        (data: IFpResult[]) => {
          return getResult(data);
        });
    return res;
  }

  async convertMolNotation(molecules: string[], targetNotation: DG.chem.Notation): Promise<string[]> {
    const t = this;
    const res =
      await this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
        t.parallelWorkers[i].convertMolNotation(segment, targetNotation),
      (data: string[][]) => {
        return ([] as string[]).concat(...data);
      });
    return res;
  }

  async getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules?: string[]):
    Promise<[RuleId, boolean[]][]> {
    const t = this;
    const fooGather = (data: {[rule in RuleId]?: boolean[]}[]): [RuleId, boolean[]][] => {
      const result: {[rule in RuleId]?: boolean[]} = {};
      for (let k = 0; k < data.length; ++k) {
        const part = data[k];
        for (const ruleId of Object.keys(part)) {
          result[ruleId as RuleId] ??= [];
          result[ruleId as RuleId] = result[ruleId as RuleId]!.concat(...part[ruleId as RuleId]!);
        }
      }
      return Object.entries(result) as [RuleId, boolean[]][];
    };
    return molecules ? this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
      t.parallelWorkers[i].getStructuralAlerts(alerts, segment), fooGather) :
      this._doParallel((i: number, _nWorkers: number) => t.parallelWorkers[i].getStructuralAlerts(alerts), fooGather);
  }

  async getRGroups(molecules: string[], coreMolecule: string, coreIsQMol: boolean, options?: string): 
    Promise<IRGroupAnalysisResult> {
    /* const t = this;
    const res = await this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
      t.parallelWorkers[i].rGroupAnalysis(segment, coreMolecule, coreIsQMol, options),
    (data: Array<IRGroupAnalysisResult>): IRGroupAnalysisResult => {
      const colNames = data[0].colNames;
      const cols = Array<Array<string>>(colNames.length).fill([]);
      for (let i = 0; i < colNames.length; i++) {
        for (let j = 0; j < data.length; j++)
          cols[i] = cols[i].concat(data[j].smiles[i]);
      }
      return {colNames: colNames, smiles: cols};
    }); */

    // R group analysis does not support parallelization, so we will use the first worker
    const res = await this.parallelWorkers[0].rGroupAnalysis(molecules, coreMolecule, coreIsQMol, options);
    return res;
  }

  async invalidateCache(): Promise<void> {
    const t = this;
    this._doParallel(
      (i: number) => {
        return t.parallelWorkers[i].invalidateCache();
      },
      (_: any) => {
        return;
      });
  }

  async setTerminateFlag(flag: boolean): Promise<void> {
    const t = this;
    this._doParallel(
      (i: number) => {
        return t.parallelWorkers[i].setTerminateFlag(flag);
      },
      (_: any) => {
        return;
      });
  }

  async mmpGetFragments(molecules: string[]): Promise<IMmpFragmentsResult> {
    const t = this;

    const getResult = (data: IMmpFragmentsResult[]): IMmpFragmentsResult => {
      return {
        frags: ([] as [string, string][][]).concat(...data.map(((it) => it.frags))),
        smiles: ([] as Array<string>).concat(...data.map(((it) => it.smiles))),
      };
    };

    const res = await this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
      t.parallelWorkers[i].mmpGetFragments(segment),
    (data: IMmpFragmentsResult[]) => {
      return getResult(data);
    });

    return res;
  }

  async mmpLinkFragments(cores: string [], fragments: string []): Promise<string[]> {
    const t = this;
    const res = await this._initParallelWorkersArray([cores, fragments], (i: number, segment: string[][]) =>
      t.parallelWorkers[i].mmpLinkFragments(segment[0], segment[1]),
    (data: string[][]): string[] => {
      return ([] as string[]).concat(...data);
    });

    return res;
  }

  async mmpGetMcs(molecules: [string, string][]): Promise<string[]> {
    const t = this;

    const res = await this._initParallelWorkers(molecules, (i: number, segment: [string, string][]) =>
      t.parallelWorkers[i].mmpGetMcs(segment),
    (data: string[][]): string[] => {
      return ([] as string[]).concat(...data);
    });
    return res;
  }

  async getInverseSubstructuresAndAlign(cores: string[], from: string[], to: string[]):
    Promise<InverseSubstructureRes> {
    const t = this;
    const res = await this._initParallelWorkersArray([cores, from, to], (i: number, segment: string[][]) =>
      t.parallelWorkers[i].getInverseSubstructuresAndAlign(segment[0], segment[1], segment[2]),
    (data: InverseSubstructureRes[]): InverseSubstructureRes => {
      return {
        inverse1: ([] as (ISubstruct | null)[]).concat(...data.map(((it) => it.inverse1))),
        inverse2: ([] as (ISubstruct | null)[]).concat(...data.map(((it) => it.inverse2))),
        fromAligned: ([] as string[]).concat(...data.map(((it) => it.fromAligned))),
        toAligned: ([] as string[]).concat(...data.map(((it) => it.toAligned))),
      };
    });
    return res;
  }

  async getMCS(molecules: string[], exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
    // MCS does not support parallelization, so we will use the first worker
    return await this.parallelWorkers[0].mostCommonStructure(molecules, exactAtomSearch, exactBondSearch);
  }
}
