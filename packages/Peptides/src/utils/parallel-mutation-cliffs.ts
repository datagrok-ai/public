import * as type from './types';
import {MutationCliffsOptions} from './algorithms';

export type ParallelMutationReturnType = {
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
    options: MutationCliffsOptions = {}): Promise<type.MutationCliffs> {
    const substitutionsInfo: type.MutationCliffs = new Map();
    try {
      const currentTargetIdx = options.targetCol?.cat!.indexOf(options.currentTarget!) ?? -1;

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
      options.targetCol?.cat && (options.targetCol.cat = options.targetCol.cat.slice());
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = Math.floor(idx * chunkSize);
          const endIdx = idx === this._workerCount - 1 ? matSize : Math.floor((idx + 1) * chunkSize);
          this._workers[idx].postMessage(
            {startIdx, endIdx, activityArray, monomerInfoArray, settings: options, currentTargetIdx});
          this._workers[idx].onmessage = ({data: {pos, seq1Idxs, seq2Idxs, error}}): void => {
            if (error) {
              this._workers[idx]?.terminate();
              rejectWorker(error);
            } else {
              this._workers[idx].terminate();
              resolveWorker({pos, seq1Idxs, seq2Idxs});
            }
          };
        });
      }

      const results = await Promise.all(promises);
      const monomerPositionsMap = new Map<string, number>();
      monomerInfoArray.forEach((monomerInfo, i) => {
        monomerPositionsMap.set(monomerInfo.name, i);
      });
      results.filter(Boolean).forEach((result) => {
        for (let i = 0; i < result.pos.length; i++) {
          //getting monomers from monomerInfoArray by position
          const monomerPos = monomerPositionsMap.get(result.pos[i])!;
          const monomer1Cat = monomerInfoArray[monomerPos].rawData[result.seq1Idxs[i]];
          const monomer1 = monomerInfoArray[monomerPos].cat![monomer1Cat];
          const monomer2Cat = monomerInfoArray[monomerPos].rawData[result.seq2Idxs[i]];
          const monomer2 = monomerInfoArray[monomerPos].cat![monomer2Cat];

          // filling map
          if (!substitutionsInfo.has(monomer1))
            substitutionsInfo.set(monomer1, new Map());
          if (!substitutionsInfo.has(monomer2))
            substitutionsInfo.set(monomer2, new Map());
          const position1Map = substitutionsInfo.get(monomer1)!;
          const position2Map = substitutionsInfo.get(monomer2)!;
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
    } catch (error) {
      this.terminate();
      console.error(error);
    }
    return substitutionsInfo;
  }

  public terminate(): void {
    this._workers?.forEach((worker) => {
      try {
        worker?.terminate();
      } catch (error) {
        console.error(error);
      }
    });
  }
}
