import {dummy} from '../cobra/dummy';
import type {CobraModelData} from '../../escher_src/src/ts/types';
import {modelFromJsonData, Solution} from './cobraSolver';

// DEPRECATED, use WorkerCobraSolver from './cobra' instead
class WorkerCobraSolver {
  private static _lastOptimizationPromise: Promise<Solution> = Promise.resolve({} as Solution);
  private static _FBAWorker: Worker | null = null;
  static async run_optimization(model_data: CobraModelData | null) {
    this._FBAWorker ??= new Worker(new URL('./cobrasolver-worker', import.meta.url));
    // make sure the model is copied
    if (!model_data)
      throw new Error('Cannot run optimization without a model loaded');
    this._lastOptimizationPromise = this._lastOptimizationPromise
      .catch(() => {
        return Promise.resolve(null);
      })
      .then(() => {
        return new Promise<Solution>((resolve, reject) => {
          try {
            this._FBAWorker!.onmessage = (ev: {data: Solution | {error: any}}) => {
              try {
                if ('error' in ev.data)
                  throw new Error(ev.data.error);
                resolve(ev.data);
              } catch (e) {
                reject(e);
              }
            };
            this._FBAWorker!.postMessage(modelFromJsonData(model_data));
          } catch (e) {
            reject(e);
          }
        });
      });
    return this._lastOptimizationPromise;
  }
}
