class RdKitParallel {
  
  constructor(nWorkers = -1) {
    if (nWorkers <= 0) {
      const cpuLogicalCores = window.navigator.hardwareConcurrency;
      this._nWorkers = Math.max(1, cpuLogicalCores - 2);
    } else {
      this._nWorkers = nWorkers;
    }
  }
  
  async init(webRoot) {
    this._workers = [];
    for (let k = 0; k < this._nWorkers; ++k) {
      let worker = new RdKitWorkerProxy(webRoot);
      await worker.moduleInit();
      this._workers.push(worker);
    }
  }
  
  async _doParallel(fooScatter, fooGather = async (d) => null) {

    let promises = [];
    const nWorkers = this._nWorkers;
    for (let i = 0; i < nWorkers; i++) {
      promises[i] = fooScatter(i, nWorkers);
    }
    let data = await Promise.all(promises);
    return fooGather(data);
    
  }
  
  async substructInit(dict) {
   
    let t = this;
    return this._doParallel(
      async (i, nWorkers) => {
        const length = dict.length;
        const segmentLength = Math.floor(length / nWorkers);
        t.segmentLength = segmentLength;
        const segment = i < (nWorkers - 1) ?
          dict.slice(i * segmentLength, (i + 1) * segmentLength) :
          dict.slice(i * segmentLength, length);
        return t._workers[i].substructInit(segment);        
      });
    
  }
  
  async substructSearch(query) {
  
    let t = this;
    return this._doParallel(
      async (i, nWorkers) => {
        return t._workers[i].substructSearch(query);        
      },
      async (data) => {
        for (let k = 0; k < data.length; ++k) {
          data[k] = JSON.parse(data[k]);
          data[k] = data[k].map(a => a + t.segmentLength * k);
        }
        return [].concat.apply([], data);
      });
  
  }

}