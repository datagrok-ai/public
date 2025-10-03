import type { CobraModelData } from '../../escher_src/src/ts/types';

export class WorkerCobraSolver {

  private static _lastOptimizationPromise: Promise<{fluxes: Float32Array, reactionNames: string[]} | null> = Promise.resolve(null);
  private static _FBAWorker: Worker | null = null;
  static async run_optimization(model_data: CobraModelData | null) {
    this._FBAWorker ??= new Worker(new URL('./glpkFBA-worker', import.meta.url));
    // make sure the model is copied
    if (!model_data)
      throw new Error('Cannot run optimization without a model loaded');
    this._lastOptimizationPromise = this._lastOptimizationPromise
      .catch(() => {
            return Promise.resolve(null);
        })
      .then(() => {
        return new Promise<{fluxes: Float32Array, reactionNames: string[]}>((resolve, reject) => {
          try {
            this._FBAWorker!.onmessage = (ev: {data: {fluxes: Float32Array, reactionNames: string[]} | {error: any}}) => {
                try {
                    if ('error' in ev.data)
                      throw new Error(ev.data.error);
                    resolve(ev.data)
                } catch (e) {
                    reject(e);
                }
            }
            this._FBAWorker!.postMessage(model_data);

          } catch (e) {
            reject(e);
          }
        })
      })
    return this._lastOptimizationPromise;
  }

  static async get_extreme_points(model_data: CobraModelData | null) {
    if (!model_data)
      throw new Error('Cannot run optimization without a model loaded');
    const reactions = model_data.reactions;
    const totalRuns = reactions.length * 2;
    const numWorkers = Math.min(Math.max((navigator.hardwareConcurrency ?? 1) - 2, 1), 20);
    const workers = new Array(numWorkers).fill(0).map(() => new Worker(new URL('./glpk-extreme-worker', import.meta.url)));
    const runsPerWorker = Math.ceil(totalRuns / numWorkers);
    const promises: Promise<{fluxes: Float32Array[], reactionNames: string[]} | {error: any}>[] = [];

    for (let w = 0; w < numWorkers; w++) {
      const start = w * runsPerWorker;
      const end = Math.min(start + runsPerWorker, totalRuns);
      if (start >= end)
        break;
      const worker = workers[w];
      const promise = new Promise<{fluxes: Float32Array[], reactionNames: string[]} | {error: any}>((resolve) => {
        worker.onmessage = (ev: {data: {fluxes: Float32Array[], reactionNames: string[]} | {error: any}}) => {
            if ('error' in ev.data)
              resolve({error: ev.data.error});
            else
              resolve({fluxes: ev.data.fluxes, reactionNames: ev.data.reactionNames});
        };
        worker.postMessage({model: model_data, start, end});
      });
      promises.push(promise);
    }
    const results = await Promise.all(promises);
    workers.forEach(worker => worker.terminate());
    const allFluxes: Float32Array[] = [];
    let reactionNames: string[] = [];
    const errors: any[] = [];
    for (const res of results) {
      if ('error' in res) {
        errors.push(res.error);
      } else {
        allFluxes.push(...res.fluxes);
        if (reactionNames.length === 0)
          reactionNames = res.reactionNames; // assuming all workers return the same reaction names
      }
    }

    if (errors.length > 0) {
      if (errors.length === results.length) {
        throw new Error(`All workers failed: ${errors.map(e => e?.toString()).join('; ')}`);
      } else {
        console.warn(`Some workers failed: ${errors.map(e => e?.toString()).join('; ')}`);
      }
    }
    return {solutions: allFluxes, reactionNames};
  }

  static async runSampling(model_data: CobraModelData | null, samplesCount: number = 1000, thinning: number = 20) {
    if (!model_data)
      throw new Error('Cannot run optimization without a model loaded');
    const worker = new Worker(new URL('./sampler-worker', import.meta.url));
    return new Promise<Float32Array>((resolve, reject) => {
      try {
        worker.onmessage = (ev: {data: Float32Array | {error: any}}) => {
            try {
                worker.terminate();
                if ('error' in ev.data)
                  throw new Error(ev.data.error);
                resolve(ev.data)
            } catch (e) {
                reject(e);
            }
        }
        worker.postMessage({model: model_data, samples: samplesCount, thinning});
      } catch (e) {
        reject(e);
      }
    })
  }

}