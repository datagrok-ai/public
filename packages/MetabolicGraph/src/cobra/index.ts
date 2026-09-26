/* eslint-disable camelcase */
import type {CobraModelData} from '../../escher_src/src/ts/types';
import type {PrecomputedWarmup} from './sampler-wrapper';
import {dummy} from './dummy';

export class WorkerCobraSolver {
  private static _lastOptimizationPromise: Promise<{fluxes: Float32Array, reactionNames: string[]} | null> = Promise.resolve(null);
  private static _FBAWorker: Worker | null = null;
  static async run_optimization(model_data: CobraModelData | null) {
    this._FBAWorker ??= new Worker(new URL('./glpk-js-fba-worker', import.meta.url));
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
                resolve(ev.data);
              } catch (e) {
                reject(e);
              }
            };
            this._FBAWorker!.postMessage(model_data);
          } catch (e) {
            reject(e);
          }
        });
      });
    return this._lastOptimizationPromise;
  }

  static async runSampling(model_data: CobraModelData | null, samplesCount: number = 1000, thinning: number = 20, precomputedWarmup?: PrecomputedWarmup) {
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
            resolve(ev.data);
          } catch (e) {
            reject(e);
          }
        };
        // a crash inside the module (e.g. an abort) never posts a message
        worker.onerror = (ev) => {
          worker.terminate();
          reject(new Error(`Sampling failed: ${ev.message}`));
        };
        worker.postMessage({model: model_data, samples: samplesCount, thinning, precomputedWarmup});
      } catch (e) {
        reject(e);
      }
    });
  }
}
