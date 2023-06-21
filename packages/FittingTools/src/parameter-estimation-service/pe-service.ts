import {PeServiceWorkerClient} from './pe-service-worker-client';

export class PeService {
  workerCount: number;
  _initWaiters: Promise<any>[];
  timesInitialized = 0;
  parallelWorkers: PeServiceWorkerClient[] = [];
  segmentLength: number = 0;

  constructor() {
    const cpuLogicalCores = window.navigator.hardwareConcurrency;
    this.workerCount = Math.max(1, cpuLogicalCores - 2);
    this._initWaiters = new Array(this.workerCount);
  }

  static async create() : Promise<PeService> {
    const service = new PeService();
    await service.init('');
    return service;
  }
  async init(webRoot: string): Promise<void> {
      let workerClient = null;
      for (let i = 0; i < this.workerCount; ++i) {
        workerClient = new PeServiceWorkerClient();
        this.parallelWorkers[i] = workerClient;
        this._initWaiters.push(workerClient.moduleInit(webRoot));
      }

    await Promise.all(this._initWaiters);
    if (this.timesInitialized++ === 0)
      console.log('RDKit Service was initialized');
  }

  async _doParallel(fooScatter: (i: number, workerCount: number) => Promise<any>, fooGather = (_: any) => []): Promise<any> {
    const promises = [];
    const workerCount = this.workerCount;
    for (let i = 0; i < workerCount; i++)
      promises[i] = fooScatter(i, workerCount);

    const data = await Promise.all(promises);
    return fooGather(data);
  }

  async fit(x: number[], y: number[], counts: number[]): Promise<string[]> {
    const t = this;
    let segmentLength = -1;
    const length = counts.length;
    return this._doParallel(
      (i: number, workerCount: number) => {
        segmentLength = Math.floor(length / workerCount);
        const segCounts = i < (workerCount - 1) ?  counts.slice(i * segmentLength, (i + 1) * segmentLength) :
                                                   counts.slice(i * segmentLength, length);
        const segXs : number[] = [];
        const segYs : number[] = [];
        let offset = 0;
        for (let k=0; k<segCounts.length; ++k) {
          for (let n = 0; n < segCounts[k]; ++n) {
            segXs.push(... x.slice(offset, offset + counts[n]));
            segYs.push(... y.slice(offset, offset + counts[n]));
            offset += counts[n];
          }
        }
       return t.parallelWorkers[i].fit(segXs, segYs, segCounts)},
      (data: any) => {
        return [].concat(...data);
      });
  }
}
