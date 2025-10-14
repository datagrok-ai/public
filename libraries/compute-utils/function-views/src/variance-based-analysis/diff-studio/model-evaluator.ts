/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DiffGrok} from '../../fitting-view';
import {OutputDataFromUI} from '../sa-outputs-routine';

import * as DGL from 'diff-grok';
import {DEFAULT_NUM, MIN_WORKERS_COUNT, WORKERS_COUNT_DOWNSHIFT,
  WorkerTask, RESULT_CODE,
  NO_ERRORS} from './defs';

const NULL = DG.FLOAT_NULL;

/** Diff Studio model analyzer */
export class ModelEvaluator {
  private diffGrok: DiffGrok;
  private funcInputs: Record<string, number>[];
  private rowIdx: number = DEFAULT_NUM;
  private argColIdx: number = DEFAULT_NUM;
  private argVal: number = DEFAULT_NUM;
  private outputNames: string[];

  constructor(
    diffGrok: DiffGrok,
    funcInputs: Record<string, number>[],
    outputsOfInterest: OutputDataFromUI) {
    this.diffGrok = diffGrok;
    this.funcInputs = funcInputs;
    this.outputNames = DGL.getOutputNames(this.diffGrok.ivp);

    if (outputsOfInterest.value.row !== null)
      this.rowIdx = outputsOfInterest.value.row;
    else {
      const name = outputsOfInterest.value.colName;
      this.argColIdx = this.outputNames.indexOf(name);

      if (this.argColIdx < 0)
        throw new Error(`Incorrect arg column name: ${name}`);

      this.argVal = outputsOfInterest.value.colValue;
    }
  }

  /** Return model evaluation */
  public async getResults(): Promise<DG.Column[]> {
    const samplesCount = this.funcInputs.length;

    // Create workers
    const nThreads = Math.min(
      Math.max(MIN_WORKERS_COUNT, navigator.hardwareConcurrency - WORKERS_COUNT_DOWNSHIFT),
      samplesCount,
    );
    const workers = new Array(nThreads)
      .fill(null)
      .map((_) => new Worker(new URL('workers/analysis.ts', import.meta.url)));

    // Create tasks for workers
    const tasks = this.getTasks(samplesCount, nThreads);

    let doneWorkers = 0;
    let percentage = 0;
    const pi = DG.TaskBarProgressIndicator.create(`Analyzing... (${percentage}%)`);

    const rawArrs: Float64Array[] = [];
    const outsCount = this.outputNames.length;

    for (let i = 0; i < outsCount; ++i)
      rawArrs.push(new Float64Array(samplesCount));

    // Run optimization
    const promises = workers.map((w, workerIdx) => {
      return new Promise<void>((resolve, reject) => {
        w.postMessage({task: tasks[workerIdx]});

        w.onmessage = (e: any) => {
          w.terminate();

          ++doneWorkers;
          percentage = Math.floor(100 * (doneWorkers + 1) / nThreads);

          if (percentage < 99.9)
            pi.update(percentage, `Analyzing... (${percentage}%)`);
          else
            pi.update(percentage, 'Analyzing... (100%)');

          if (e.data.callResult === RESULT_CODE.SUCCEED) {
            const rowVals = e.data.rowVals;

            e.data.analysisRes.forEach((res: string, rowIdx: number) => {
              if (res === NO_ERRORS) {
                for (let k = 0; k < outsCount; ++k)
                  rawArrs[k][tasks[workerIdx].indices[rowIdx]] = rowVals[rowIdx][k];
              } else {
                for (let k = 0; k < outsCount; ++k)
                  rawArrs[k][tasks[workerIdx].indices[rowIdx]] = NULL;
              }
            });
          } else {
            reject(e.data.msg ?? 'error in webworker sensitivity analysis');
            return;
          }
          resolve();
        };

        w.onerror = (e) => {
          w.terminate();
          console.error(e);
          reject(e);
        };
      });
    }); // promises

    await Promise.all(promises);

    pi.close();

    return this.outputNames.map((name, idx) => DG.Column.fromFloat64Array(name, rawArrs[idx]));
  }

  /** Return  tasks splitted into batches */
  private getTasks(samplesCount: number, nThreads: number): WorkerTask[] {
    const tasks: WorkerTask[] = [];
    let inputIdx = 0;
    let inputVec: Float64Array;
    const chunkSize = Math.floor(samplesCount / nThreads);
    let remainder = samplesCount % nThreads;
    let batchSize: number;
    let extra: number;

    for (let taskIdx = 0; taskIdx < nThreads; ++taskIdx) {
      const pipelines: DGL.Pipeline[] = [];
      const inputVecs: Float64Array[] = [];
      const indices: number[] = [];
      extra = remainder > 0 ? 1 : 0;
      batchSize = chunkSize + extra;
      --remainder;

      for (let cur = 0; cur < batchSize; ++cur) {
        inputVec = DGL.getInputVector(this.funcInputs[inputIdx], this.diffGrok.ivp);
        inputVecs.push(inputVec);
        pipelines.push(this.diffGrok.pipelineCreator.getPipeline(inputVec));
        indices.push(inputIdx);
        ++inputIdx;

        if (inputIdx >= samplesCount)
          break;
      }

      tasks.push({
        ivp2ww: this.diffGrok.ivpWW,
        pipelines: pipelines,
        indices: indices,
        inputVecs: inputVecs,
        rowIdx: this.rowIdx,
        argColIdx: this.argColIdx,
        argVal: this.argVal,
      });
    }

    return tasks;
  } // getTasks
}; // ModelEvaluator
