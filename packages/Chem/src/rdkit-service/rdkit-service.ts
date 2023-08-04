import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import {Fingerprint} from '../utils/chem-common';
import {RuleId} from '../panels/structural-alerts';
import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {IFpResult} from './rdkit-service-worker-similarity';
import {LockedEntity} from '../utils/locked-entitie';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {getRdKitModule} from '../package';
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
};

export class RdKitService {
  workerCount: number;
  _initWaiters?: Promise<any>[];
  timesInitialized = 0;
  parallelWorkers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;
  moleculesSegmentsLengths: Uint32Array;
  segmentLengthPatternFp: number = 0;
  moleculesSegmentsLengthsPatternFp: Uint32Array;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this.workerCount = Math.max(1, cpuLogicalCores - 2);
    this.moleculesSegmentsLengths = new Uint32Array(this.workerCount);
    this.moleculesSegmentsLengthsPatternFp = new Uint32Array(this.workerCount);
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
   * @param {function (workerIdx: number, workerCount: number): Promise<TMap>} map - splits the data by number of workers
   * and calls function from rdkit service worker client (basicaly action which we need to perform inside worker -
   * getFingerprints, searchSubstructure etc.)
   * @param {function (_: TMap[]): TReduce} reduce - fucntion which combines results collected from web workers into single result
   * */
  async _doParallel<TMap, TReduce>(
    map: (workerIdx: number, workerCount: number) => Promise<TMap>,
    reduce: (_: TMap[]) => TReduce = (_: TMap[]) => [] as TReduce): Promise<TReduce> { //(_: TMap[]) => [] as TReduce is a default function
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
    pogressFunc: (progress: number) => void): Promise<IParallelBatchesRes> {
      let terminateFlag = false;

      const setTerminateFlag = () => { terminateFlag = true; }
      const getTerminateFlag = () => { return terminateFlag; }
      const t = this;
      const dataLength = data.length;
      const increment = 50;
      const workingIndexes = new Array<{start: number, end: number, increment: number}>(this.workerCount).fill({start: 0, end: 0, increment});
      const sizePerWorker = Math.floor(dataLength / this.workerCount);
      
      // distribute data between workers, makes sure that same molecules always end up in the same worker, important for caching
      for (let i = 0; i < this.workerCount; i++) {
        workingIndexes[i] = {start: i * sizePerWorker, end: i === this.workerCount - 1 ? dataLength : (i + 1) * sizePerWorker, increment};
      }
      const lockedCounter = new LockedEntity(0);
      
      let processedMolecules = 0;
      let moleculesPerProgress = Math.min(Math.max(Math.floor(dataLength / 100), 10), 100);
      let nextProgressCheck = moleculesPerProgress;
      const promises = t.parallelWorkers.map((_, idx) => {
          const post = async () => {
              if (workingIndexes[idx].start >= workingIndexes[idx].end || terminateFlag) {
                  return;
              }
              const part = data.slice(workingIndexes[idx].start, Math.min(workingIndexes[idx].end, workingIndexes[idx].start + workingIndexes[idx].increment));
              const batchResult = await map(part, idx, t.parallelWorkers.length, workingIndexes[idx].start);
              updateRes(batchResult, res, part.length, workingIndexes[idx].start);
              workingIndexes[idx].start += part.length;
              
              await lockedCounter.unlockPromise();
              lockedCounter.lock();
              processedMolecules = lockedCounter.value;
              const end = Math.min(processedMolecules + workingIndexes[idx].increment, dataLength);
              if (processedMolecules >= nextProgressCheck) {
                nextProgressCheck += moleculesPerProgress;
                moleculesPerProgress *= 1.5;
                // increment *= 1.2;
                //increment = Math.floor(increment);
                pogressFunc(processedMolecules/dataLength);
              }
              workingIndexes[idx].increment = Math.floor(workingIndexes[idx].increment * 1.05);
              lockedCounter.value = end;
              lockedCounter.release();
              await post();
          }      
          return post();    
  })
  const getProgress = () => processedMolecules / dataLength;
  return {getProgress, setTerminateFlag, getTerminateFlag, promises};      
}


   /**
   * Calls _doParallel with pre-defined map function which splits data by number of workers
   * @async
   * @param {string[]} molecules - list of molecules to split by workers
   * @param {function (workerIdx: number, workerCount: number): Promise<TMap>} workerFunc - function from rdkit service worker
   *  client (basicaly action which we need to perform inside worker - getFingerprints, searchSubstructure etc.)
   * @param {function (_: TMap[]): TReduce} reduce - function which combines results collected from web workers into single result
   * */
  async _initParallelWorkers<TMap, TReduce>(molecules: string[],
    workerFunc: (workerIdx: number, moleculesSegment: string[]) => Promise<TMap>,
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
   * @param {string} queryMolBlockFailover - smart to filter by (is used if creation of RDMol object from query parameter failed)
   * @param {SubstructureSearchWithFpResult} result - object where to put the results of substructure search and created fp
   * @param {function (progress: number): void} progressFunc - function which should be run on each progress iteration
   * @param {boolean} molecules - list of molecules to search
   * @param {boolean} createSmiles - whether canonical smiles should be created besides fp (should be true if molecules are not canonical smiles)
   * */
  async searchSubstructureWithFps(query: string, queryMolBlockFailover: string, result: SubstructureSearchWithFpResult,
    progressFunc: (progress: number) => void, molecules: string[], createSmiles = false) {
      const queryMol = getQueryMolSafe(query, queryMolBlockFailover, getRdKitModule());
      if (!queryMol)
        throw new Error(`Chem | Invalid search pattern: ${query}`);
      let fpRdKit: Uint8Array;
      try {
        fpRdKit = queryMol.get_pattern_fp_as_uint8array();
      } catch (e: any) {
        throw new Error(`Chem | Substructure Search failed with error: ${e.toString()}`);
      }
      

      const updateRes = (batchRes: SubstructureSearchBatchResult, res: SubstructureSearchWithFpResult, length: number, index: number) => {
        if (!res.fpsRes) {
          res.fpsRes = {
          fps: new Array<Uint8Array | null>(molecules.length).fill(null),
          smiles: createSmiles ? new Array<string | null>(molecules.length).fill(null) : null
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
        }

        return this._doParallelBatches(molecules, result, updateRes, async (batch, workerIdx, _workerCount, batchStartIdx) => {
          let fpResult: IFpResult;
          if (!result.fpsRes || !result.fpsRes.fps[batchStartIdx])
            fpResult = await this.parallelWorkers[workerIdx].getFingerprints(Fingerprint.Pattern, batch, createSmiles);
          else {
            fpResult = {
              fps: result.fpsRes.fps.slice(batchStartIdx, batchStartIdx + batch.length),
              smiles: createSmiles ? result.fpsRes.smiles!.slice(batchStartIdx, batchStartIdx + batch.length) : batch
            }
          }

          if(query === '')
            return {matches: new BitArray(batch.length, true), fpRes: fpResult};
          // *********** FILTERING using fingerprints
          const patternFpUint8Length = 256;
          const patternFpFilterBitArray = new BitArray(batch.length, false);
            try {
              
              checkEl:
              for (let i = 0; i < batch.length; ++i) {
                if (fpResult.fps[i]) {
                  for (let j = 0; j < patternFpUint8Length; ++j) {
                    if ((fpResult.fps[i]![j] & fpRdKit[j]) != fpRdKit[j])
                      continue checkEl;
                  }
                  patternFpFilterBitArray.setFast(i, true);
                }
              }
            } catch (e: any) { 
             throw new Error(`Chem | Substructure Search failed with error: ${e.toString()}`);
            }

            const filteredMolecules = Array<string>(patternFpFilterBitArray.trueCount());
            let counter = 0;
            for (let i = -1; (i = patternFpFilterBitArray.findNext(i)) !== -1;) {
              filteredMolecules[counter] = createSmiles ? fpResult.smiles![i]! : batch[i];
              counter++;
            }

          // *********** DONE FILTERING using fingerprints
          // filter using substruct search on already prefiltered dataset
          const substructRes: Uint32Array = await this.parallelWorkers[workerIdx]
            .searchSubstructure(query, queryMolBlockFailover, filteredMolecules);
          
          const matchesBitArray = BitArray.fromUint32Array(filteredMolecules.length, substructRes);
          // restore the indexes of prefiltered molecules on the whole dataset
          const restoredBitArray = new BitArray(batch.length, false);
          let matchesCounter = 0;
          for (let i = -1; (i = patternFpFilterBitArray.findNext(i)) != -1;) {
            if (matchesBitArray.getBit(matchesCounter))
            restoredBitArray.setBit(i, true);
            matchesCounter++;
          }
          return {
            matches: restoredBitArray,
            fpRes: fpResult
          }
        }, progressFunc)
    }




   /**
   * Returns fingerprints for provied molecules
   * @async
   * @param {Fingerprint} fingerprintType - type of fingerprint (Morgan or Pattern)
   * @param {boolean} molecules - optional parameter, if passed RDMols objects are created on the fly (array of
   * predefined RDMols is not used) 
   * @param {boolean} getCanonicalSmiles - optional parameter, if passed, in addition to fps function returns canonical smiles
   * */
  async getFingerprints(fingerprintType: Fingerprint, molecules?: string[], getCanonicalSmiles?: boolean): Promise<IFpResult> {
    const t = this;
    const getResult = (data: IFpResult[]): IFpResult => {
      return {
        fps: ([] as Array<Uint8Array | null>).concat(...data.map((it => it.fps))),
        smiles: getCanonicalSmiles ? ([] as Array<string | null>).concat(...data.map((it => it.smiles))) : null
      };
    }
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


  async convertMolNotation(targetNotation: string): Promise<string[]> {
    const t = this;
    return this._doParallel(
      (i: number, nWorkers: number) => {
        return t.parallelWorkers[i].convertMolNotation(targetNotation);
      },
      (data: any) => {
        for (let k = 0; k < data.length; ++k)
          data[k] = data[k].map((a: number) => a + t.segmentLength * k);

        return [].concat(...data);
      });
  }

  async getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules?: string[]): Promise<[RuleId, boolean[]][]> {
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
}
