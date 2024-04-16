/** Monte Carlo optimizer */

import {Extremum} from './optimizer-misc';
import {optimizeNM} from './optimizer-nelder-mead';
import {sampleParams} from './optimizer-sampler';

export async function optimize(
  objectiveFunc: (x: Float32Array) => Promise<number>,
  paramsBottom: Float32Array,
  paramsTop: Float32Array,
  samplesCount: number): Promise<Extremum[]> {
  const params = sampleParams(samplesCount, paramsTop, paramsBottom);

  // const threadCount = Math.max(navigator.hardwareConcurrency - 2, 1);

  // const workers = new Array(threadCount).fill(null)
  //   .map(() => new Worker(new URL('./optimizer-worker', import.meta.url)));
  // const chunkSize = samplesCount / threadCount;
  // const promises = new Array<Promise<{extremums: Extremum[], warnings: string[]}>>(threadCount);

  // for (let i = 0; i < threadCount; i++) {
  //   const start = Math.floor(i * chunkSize);
  //   const end = (i === threadCount - 1) ? samplesCount : Math.floor((i + 1) * chunkSize);
  //   workers[i].postMessage({objectiveFunc, params, start, end});
  //   promises[i] = new Promise<{extremums: Extremum[], warnings: string[]}>((resolveWorker) => {
  //     workers[i].onmessage = ({data: {extremums, warnings}}): void => {
  //       resolveWorker({extremums, warnings});
  //     };
  //   });
  // }
  // const resArray = await Promise.all(promises);

  // let extremumsData: Extremum[] = [];
  // let warningsAll: string [] = [];

  // resArray.forEach((resItem) => {
  //   extremumsData = extremumsData.concat(...resItem.extremums);
  //   warningsAll = warningsAll.concat(...resItem.warnings);
  // });

  // setTimeout(() => {
  //   workers.forEach((worker) => {
  //     worker.terminate();
  //   });
  // }, 0);

  const extremums = new Array<Extremum>(samplesCount);
  for (let i = 0; i < samplesCount; i++)
    extremums[i] = await optimizeNM(objectiveFunc, params[i]);

  console.log(extremums);

  return extremums;
}
