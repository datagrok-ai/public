import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Fingerprint} from '../utils/chem-common';
import { RuleId } from '../panels/structural-alerts';
import * as DG from 'datagrok-api/dg';

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

  async _doParallel(
    map: (i: number, workerCount: number) => Promise<any>,
    reduce: (_: any) => any = (_: any) => []): Promise<any> {
    const promises = [];
    const workerCount = this.workerCount;
    for (let i = 0; i < workerCount; i++)
      promises[i] = map(i, workerCount);

    const data = await Promise.all(promises);
    return reduce(data);
  }

  async _initParallelWorkers(molecules: string[], func: any, postFunc: any): Promise<any> {
    const t = this;
    return this._doParallel(
      (i: number, nWorkers: number) => {
        const length = molecules.length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segment = i < (nWorkers - 1) ?
          molecules.slice(i * segmentLength, (i + 1) * segmentLength) :
          molecules.slice(i * segmentLength, length);
        t.moleculesSegmentsLengths![i] = segment.length;
        return func(i, segment);
      },
      postFunc,
    );
  }

  async initMoleculesStructures(molecules: string[], useSubstructLib?: boolean)
    : Promise<number> {
    return this._initParallelWorkers(molecules, (i: number, segment: any) =>
      this.parallelWorkers[i].initMoleculesStructures(segment, useSubstructLib),
    () => {});
  }

  async searchSubstructure(query: string, queryMolBlockFailover: string, molecules?: string[], 
    useSubstructLib?: boolean): Promise<DG.BitSet> {
    const t = this;
    if (molecules) { //need to recalculate segments lengths since when using pattern fp numbsr of molecules to search can be less that totla molecules in column
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
        const bitArray = DG.BitSet.create(totalLength);
        let counter = 0;
        for (let i = 0; i < segmentsLengths.length; i++) {
          for (let j = 0; j < segmentsLengths[i]; j++) {
            bitArray.set(counter++, !!(data[i][0] >> j & 1));
          }
        }
        return bitArray;
      });
  }


  async getFingerprints(fingerprintType: Fingerprint, molecules?: string[]): Promise<Uint8Array[]> {
    const t = this;
    const res = molecules ?
      await this._initParallelWorkers(molecules, (i: number, segment: string[]) =>
        t.parallelWorkers[i].getFingerprints(fingerprintType, segment),
      (data: Array<Uint8Array | null>[][]) => {
        return ([] as Array<Uint8Array | null>[]).concat(...data);
      }) :
      await this._doParallel(
        (i: number, _: number) => {
          return t.parallelWorkers[i].getFingerprints(fingerprintType, molecules);
        },
        (data: Array<Uint8Array | null>[][]) => {
          return ([] as Array<Uint8Array | null>[]).concat(...data);
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
