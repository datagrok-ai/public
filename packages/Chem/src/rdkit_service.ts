import {RdKitServiceWorkerClient} from './rdkit_service_worker_client';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

export class RdKitService {
  readonly _nJobWorkers = 1; // only 1 for now
  _nParallelWorkers: number;
  _jobWorkers: RdKitServiceWorkerClient[] = [];
  _parallelWorkers: RdKitServiceWorkerClient[] = [];
  _jobWorker: RdKitServiceWorkerClient | undefined;
  segmentLength: number = 0;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this._nParallelWorkers = Math.max(1, cpuLogicalCores - 2);
  }

  async init(webRoot: string): Promise<void> {
    this._parallelWorkers = [];
    this._jobWorkers = [];
    const initWaiters = [];
    for (let i = 0; i < this._nParallelWorkers; ++i) {
      const workerClient = new RdKitServiceWorkerClient();
      if (i < this._nJobWorkers) {
        this._jobWorkers[i] = workerClient;
      }
      this._parallelWorkers[i] = workerClient;
      initWaiters.push(workerClient.moduleInit(webRoot));
    }
    await Promise.all(initWaiters);
    this._jobWorker = this._jobWorkers[0];
  }

  async _doParallel(fooScatter: any, fooGather = (d: any) => []): Promise<any> {
    const promises = [];
    const nWorkers = this._nParallelWorkers;
    for (let i = 0; i < nWorkers; i++) {
      promises[i] = fooScatter(i, nWorkers);
    }
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

  async initMoleculesStructures(dict: string[]): Promise<any> {
    return this._initParallelWorkers(dict, (i: number, segment: any) =>
      this._parallelWorkers[i].initMoleculesStructures(segment),
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
