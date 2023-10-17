
import * as type from './types';

export type ParallelMutationReturnType = {
    monomers1: string[],
    monomers2: string[],
    pos: string[],
    seq1Idxs: Uint32Array,
    seq2Idxs: Uint32Array,
}
export class ParallelMutationCliffs {
  private _workers: Worker[];
  private _workerCount: number;

  constructor() {
    const threadCount = navigator.hardwareConcurrency;
    this._workerCount = Math.max(threadCount - 2, 1);
    this._workers = new Array(this._workerCount).fill(null)
      .map(() => new Worker(new URL('../workers/mutation-cliffs-worker', import.meta.url)));
  }

  public async calc(activityArray: type.RawData, monomerInfoArray: type.RawColumn[],
    settings: type.PeptidesSettings = {},
    targetOptions: {targetCol?: type.RawColumn | null, currentTarget?: string | null} = {},
  ): Promise<type.MutationCliffs> {
    console.time('parallel');
    const currentTargetIdx = targetOptions.targetCol?.cat!.indexOf(targetOptions.currentTarget!) ?? -1;

    const len = activityArray.length;
    const promises = new Array<Promise<ParallelMutationReturnType>>(this._workerCount);
    const matSize = len * (len - 1) / 2; // size of reduced upper triangular matrix
    this._workerCount = Math.min(this._workerCount, matSize);
    const chunkSize = matSize / this._workerCount;
    // monomerInfoArray[m].cat and targetCol can contain some function from datagrok-api,
    //which the worker can't serialize and fails. so we need to remove it
    monomerInfoArray.forEach((monomerInfo) => {
      monomerInfo.cat = monomerInfo.cat?.slice();
    });
    targetOptions.targetCol?.cat && (targetOptions.targetCol.cat = targetOptions.targetCol.cat.slice());
    for (let idx = 0; idx < this._workerCount; idx++) {
      promises[idx] = new Promise((resolveWorker, rejectWorker) => {
        const startIdx = Math.floor(idx * chunkSize);
        const endIdx = idx === this._workerCount - 1 ? matSize : Math.floor((idx + 1) * chunkSize);
        this._workers[idx].postMessage(
          {startIdx, endIdx, activityArray, monomerInfoArray, settings, currentTargetIdx, targetOptions});
        this._workers[idx].onmessage = ({data: {monomers1, monomers2, pos, seq1Idxs, seq2Idxs, error}}): void => {
          if (error) {
            this._workers[idx]?.terminate();
            rejectWorker(error);
          } else {
            this._workers[idx].terminate();
            resolveWorker({monomers1, monomers2, pos, seq1Idxs, seq2Idxs});
          }
        };
      });
    }

    const results = await Promise.all(promises);
    console.timeEnd('parallel');

    const substitutionsInfo: type.MutationCliffs = new Map();
    console.time('fill map');
    results.forEach((result) => {
      for (let i = 0; i< result.pos.length; i++) {
        if (!substitutionsInfo.has(result.monomers1[i]))
          substitutionsInfo.set(result.monomers1[i], new Map());
        if (!substitutionsInfo.has(result.monomers2[i]))
          substitutionsInfo.set(result.monomers2[i], new Map());
        const position1Map = substitutionsInfo.get(result.monomers1[i])!;
        const position2Map = substitutionsInfo.get(result.monomers2[i])!;
        if (!position1Map.has(result.pos[i]))
          position1Map.set(result.pos[i], new Map());
        if (!position2Map.has(result.pos[i]))
          position2Map.set(result.pos[i], new Map());
        const indexes1Map = position1Map.get(result.pos[i])!;
        const indexes2Map = position2Map.get(result.pos[i])!;
        if (!indexes1Map.has(result.seq1Idxs[i]))
          indexes1Map.set(result.seq1Idxs[i], []);
        if (!indexes2Map.has(result.seq2Idxs[i]))
          indexes2Map.set(result.seq2Idxs[i], []);
        const indexes1 = indexes1Map.get(result.seq1Idxs[i])!;
        const indexes2 = indexes2Map.get(result.seq2Idxs[i])!;
        (indexes1 as number[]).push(result.seq2Idxs[i]);
        (indexes2 as number[]).push(result.seq1Idxs[i]);
      }
    });
    console.timeEnd('fill map');
    return substitutionsInfo;
  }

  public terminate(): void {
    this._workers?.forEach((worker) => worker?.terminate());
  }
}
