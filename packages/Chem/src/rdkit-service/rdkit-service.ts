import {RdKitServiceWorkerClient} from './rdkit-service-worker-client';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Fingerprint} from '../utils/chem-common';
import {BitSet} from 'datagrok-api/dg';

export class RdKitService {
  workerCount: number;
  _initWaiters?: Promise<any>[];
  timesInitialized = 0;
  parallelWorkers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this.workerCount = Math.max(1, cpuLogicalCores - 2);
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
    fooScatter: (i: number, workerCount: number) => Promise<any>,
    fooGather = (_: any) => []): Promise<any> {
    const promises = [];
    const workerCount = this.workerCount;
    for (let i = 0; i < workerCount; i++)
      promises[i] = fooScatter(i, workerCount);

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

  async initMoleculesStructures(dict: string[])
    : Promise<any> {
    return this._initParallelWorkers(dict, (i: number, segment: any) =>
      this.parallelWorkers[i].initMoleculesStructures(segment),
    () => {});
  }

  async searchSubstructure(query: string, queryMolBlockFailover: string, bitset?: BitArray): Promise<number[]> {
    const t = this;
    return this._doParallel(
      (i: number, nWorkers: number) => {
        return bitset ?
          t.parallelWorkers[i].searchSubstructure(query, queryMolBlockFailover, i < (nWorkers - 1) ?
            bitset.getRangeAsList(i * this.segmentLength, (i + 1) * this.segmentLength) :
            bitset.getRangeAsList(i * this.segmentLength, bitset.length)) :
          t.parallelWorkers[i].searchSubstructure(query, queryMolBlockFailover);
      },
      (data: any) => {
        for (let k = 0; k < data.length; ++k) {
          // check for parse neccessity
          data[k] = JSON.parse(data[k]);
          data[k] = data[k].map((a: number) => a + t.segmentLength * k);
        }
        return [].concat(...data);
      });
  }

  async getFingerprints(fingerprintType: Fingerprint): Promise<Uint8Array[]> {
    const t = this;
    return (await this._doParallel(
      (i: number, _: number) => {
        return t.parallelWorkers[i].getFingerprints(fingerprintType);
      },
      (data: any) => {
        return [].concat(...data);
      })).map(
      (obj: any) =>
      // We deliberately choose Uint32Array over DG.BitSet here
        obj.data);
  }
}
