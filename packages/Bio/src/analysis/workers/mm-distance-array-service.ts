import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

export async function calculateMMDistancesArray(
  macromoleculeCol: DG.Column, templateIdx: number
): Promise<Float32Array> {
  const values = macromoleculeCol.toList();
  if (macromoleculeCol.semType !== DG.SEMTYPE.MACROMOLECULE)
    throw new Error('Column has to be of macromolecule type');
  const uh = UnitsHandler.getOrCreate(macromoleculeCol);
  const fnName = uh.getDistanceFunctionName();
  const threadCount = Math.min(Math.max(navigator.hardwareConcurrency - 2, 1), values.length);
  const workers = new Array(threadCount).fill(null).map((_i) =>
    new Worker(new URL('mm-distance-array-worker', import.meta.url)));
  const res = new Float32Array(values.length);
  let lmin = 0;
  let lmax = Number.MIN_VALUE;
  const promises = workers.map((worker, i) => {
    const start = Math.floor(i * values.length / threadCount);
    const end = i === workers.length - 1 ? Math.floor((i + 1) * values.length / threadCount) : values.length;
    return new Promise<void>((resolve, reject) => {
      worker.onmessage = ({data: {error, distanceArrayData, min, max}}) => {
        if (error) {
          reject(error);
        } else {
          lmin = Math.min(lmin, min);
          lmax = Math.max(lmax, max);
          res.set(distanceArrayData, start);
          resolve();
        }
      };
      worker.postMessage({fnName, values, templateIdx, start, end});
    });
  });

  try {
    await Promise.all(promises);
    res.forEach((value, index) => { res[index] = (value - lmin) / (lmax - lmin); });
    workers.forEach((worker) => worker.terminate());
  } catch (e) {
    workers.forEach((worker) => worker.terminate());
    throw e;
  }

  return res;
}
