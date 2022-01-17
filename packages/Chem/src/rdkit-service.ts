import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {chemBeginCriticalSection, chemEndCriticalSection} from "./chem-common";

export class RdKitService {
  _nParallelWorkers: number;
  _initWaiters?: Promise<any>[];
  _timesInitialized = 0;
  _parallelWorkers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this._nParallelWorkers = Math.max(1, cpuLogicalCores - 2);
  }

  async init(webRoot: string) {
    if (!this._initWaiters) {
      this._initWaiters = [];
      for (let i = 0; i < this._nParallelWorkers; ++i) {
        const workerClient = new RdKitServiceWorkerClient();
        this._parallelWorkers[i] = workerClient;
        this._initWaiters.push(workerClient.moduleInit(webRoot));
      }
    }
    await Promise.all(this._initWaiters);
    if (this._timesInitialized++ === 0)
      console.log('RDKit Service was initialized');
  }

  async _doParallel(fooScatter: any, fooGather = (d: any) => []): Promise<any> {
    const promises = [];
    const nWorkers = this._nParallelWorkers;
    for (let i = 0; i < nWorkers; i++)
      promises[i] = fooScatter(i, nWorkers);

    const data = await Promise.all(promises);
    return fooGather(data);
  }

  async _initParallelWorkers(dict: string[], func: any, postFunc: any): Promise<any> {
    const t = this;
    return this._doParallel(
      (i: number, nWorkers: number) => {
        const length = dict.length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segment = i < (nWorkers - 1) ?
          dict.slice(i * segmentLength, (i + 1) * segmentLength) :
          dict.slice(i * segmentLength, length);
        return func(i, segment);
      },
      postFunc,
    );
  }

  async initMoleculesStructures(dict: string[], usePatternFingerprints: boolean = false): Promise<any> {
    return this._initParallelWorkers(dict, (i: number, segment: any) =>
      this._parallelWorkers[i].initMoleculesStructures(segment, usePatternFingerprints),
    (resultArray: any[]) => resultArray.reduce((acc: any, item: any) => {
      item = item || {molIdxToHash: [], hashToMolblock: {}};
      return {
        molIdxToHash: [...acc.molIdxToHash, ...item.molIdxToHash],
        hashToMolblock: {...acc.hashToMolblock, ...item.hashToMolblock},
      };
    }, {molIdxToHash: [], hashToMolblock: {}}));
  }

  async searchSubstructure(query: string, querySmarts: string) {
    const t = this;
    return this._doParallel(
      (i: number, nWorkers: number) => {
        return t._parallelWorkers[i].searchSubstructure(query, querySmarts);
      },
      (data: any) => {
        for (let k = 0; k < data.length; ++k) {
          data[k] = JSON.parse(data[k]);
          data[k] = data[k].map((a: number) => a + t.segmentLength * k);
        }
        return [].concat.apply([], data);
      });
  }

  async initMorganFingerprints() {
    return this._doParallel(
      (i: number, nWorkers: number) => {
        return this._parallelWorkers[i].initMorganFingerprints();
      }, (_: any) => {
        return [];
      },
    );
  }

  async getMorganFingerprints() {
    const t = this;
    return (await this._doParallel(
      (i: number, nWorkers: number) => {
        return t._parallelWorkers[i].getMorganFingerprints();
      },
      (data: any) => {
        return [].concat.apply([], data);
      })).map(
      (obj: any) =>
      // We deliberately choose Uint32Array over DG.BitSet here
      //@ts-ignore
        new BitArray(new Uint32Array(obj.data), obj.length));
  }
}
