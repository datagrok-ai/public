import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import {Fingerprint} from '../utils/chem-common';
import { RuleId } from '../panels/structural-alerts';
import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { IFpResult } from './rdkit-service-worker-similarity';
import * as grok from 'datagrok-api/grok';
import { TERMINATE_CURRENT_SEARCH } from '../constants';

export class RdKitService {
  workerCount: number;
  _initWaiters?: Promise<any>[];
  timesInitialized = 0;
  parallelWorkers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;
  moleculesSegmentsLengths: Uint32Array;
  segmentLengthPatternFp: number = 0;
  moleculesSegmentsLengthsPatternFp: Uint32Array;
  terminateFlag: boolean = false;

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
   * Fills array of mols or SubstructLibrary object in each worker
   * @async
   * @param {string[]} molecules - list of molecules to save in each worker
   * @param {boolean} useSubstructLib - optional parameter, if set to true SubstructLibrary object if filled with mols,
   * otherwise array of mols is used
   * */
  async initMoleculesStructures<TReduce>(molecules: string[], useSubstructLib?: boolean)
    : Promise<TReduce | void> {
    return this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
      this.parallelWorkers[i].initMoleculesStructures(segment, useSubstructLib),
    () => {});
  }

   /**
   * Filters molecules by substructure
   * @async
   * @param {string} query - smiles/molblock to filter by
   * @param {string} queryMolBlockFailover - smart to filter by (is used if creation of RDMol object from query parameter failed)
   * @param {boolean} molecules - optional parameter, if passed RDMols objects are created on the fly (neither array of
   * predefined RDMols nor SubstructLibrary are used) 
   * @param {boolean} useSubstructLib - optional parameter, if set to true SubstructLibrary is used for search
   * */
  async searchSubstructure(query: string, queryMolBlockFailover: string, molecules?: string[], 
    useSubstructLib?: boolean): Promise<BitArray> {
    const t = this;

    grok.events.onCustomEvent(TERMINATE_CURRENT_SEARCH).subscribe(() => {
      this.parallelWorkers.forEach((worker) => worker.postTerminationFlag(true));
    });

    if (molecules) { //need to recalculate segments lengths since when using pattern fp number of molecules to search can be less that totla molecules in column
      this.segmentLengthPatternFp = Math.floor(molecules.length / this.workerCount);
      for (let j = 0; j < this.workerCount; j++) {
        this.moleculesSegmentsLengthsPatternFp[j] = j < (this.workerCount - 1) ? this.segmentLengthPatternFp :
          molecules.length - this.segmentLengthPatternFp * j;
      }
    }
    return this._doParallel(
      (i: number, nWorkers: number) => {
        return molecules ?
          t.parallelWorkers[i].searchSubstructure(query, queryMolBlockFailover, i < (nWorkers - 1) ?
            molecules.slice(i * this.segmentLengthPatternFp, (i + 1) * this.segmentLengthPatternFp) :
            molecules.slice(i * this.segmentLengthPatternFp, molecules.length)) :
          t.parallelWorkers[i].searchSubstructure(query, queryMolBlockFailover, undefined, useSubstructLib);
      },
      (data: Array<Uint32Array>) => {
        const segmentsLengths = molecules ? t.moleculesSegmentsLengthsPatternFp : t.moleculesSegmentsLengths;
        const totalLength = segmentsLengths.reduce((acc, cur) => acc + cur, 0);
        const bitset = new BitArray(totalLength);
        let counter = 0;
        for (let i = 0; i < segmentsLengths.length; i++) {
          for (let j = 0; j < segmentsLengths[i]; j++) {
            bitset.setBit(counter++, !!(data[i][Math.floor(j / 32)] >> j % 32 & 1));
          }
        }
        return bitset;
      });
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
}
