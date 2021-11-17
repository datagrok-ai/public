import {RdKitServiceWorkerClient} from './rdkit_service_worker_client';

export class RdKitService {

  _nWorkers: number;
  _workers: RdKitServiceWorkerClient[] = [];
  segmentLength: number = 0;

  constructor(nWorkers = -1) {
    if (nWorkers <= 0) {
      const cpuLogicalCores = window.navigator.hardwareConcurrency;
      this._nWorkers = Math.max(1, cpuLogicalCores - 2);
    } else {
      this._nWorkers = nWorkers;
    }
  }
  
  async init(webRoot: string): Promise<void> {
    this._workers = [];
    let initWaiters = [];
    for (let k = 0; k < this._nWorkers; ++k) {
      let worker = new RdKitServiceWorkerClient();
      initWaiters.push(worker.moduleInit(webRoot));
      this._workers.push(worker);
    }
    await Promise.all(initWaiters);
  }
  
  async _doParallel(fooScatter: any, fooGather = async (d: any) => []): Promise<any> {

    let promises = [];
    const nWorkers = this._nWorkers;
    for (let i = 0; i < nWorkers; i++) {
      promises[i] = fooScatter(i, nWorkers);
    }
    let data = await Promise.all(promises);
    return fooGather(data);
    
  }
  
  async substructInit(dict: string[]): Promise<any> {
   
    let t = this;
    return this._doParallel(
      async (i: number, nWorkers: number) => {
        const length = dict.length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segment = i < (nWorkers - 1) ?
          dict.slice(i * segmentLength, (i + 1) * segmentLength) :
          dict.slice(i * segmentLength, length);
        return t._workers[i].substructInit(segment);
      },
      async (resultArray) => resultArray.reduce((acc: any, item: any) => {
        item = item || { molIdxToHash: [], hashToMolblock: {} };
        return {
          molIdxToHash: [ ...acc.molIdxToHash, ...item.molIdxToHash ],
          hashToMolblock: { ...acc.hashToMolblock, ...item.hashToMolblock }
        }
      }, { molIdxToHash: [], hashToMolblock: {} })
    );
  }
  
  async substructSearch(query: string, querySmarts: string) {
  
    let t = this;
    return this._doParallel(
      async (i: number, nWorkers: number) => {
        return t._workers[i].substructSearch(query, querySmarts);        
      },
      async (data: any) => {
        for (let k = 0; k < data.length; ++k) {
          data[k] = JSON.parse(data[k]);
          data[k] = data[k].map((a: number) => a + t.segmentLength * k);
        }
        return [].concat.apply([], data);
      });
  
  }

}