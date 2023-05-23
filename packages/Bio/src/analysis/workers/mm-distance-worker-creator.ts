import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

export async function calcMmDistanceMatrix(column: DG.Column<any>): Promise<Float32Array> {
  const values = column.toList();
  const worker = new Worker(new URL('./mm-distance-worker.ts', import.meta.url));
  if (column.semType !== DG.SEMTYPE.MACROMOLECULE)
    throw new Error('Column has to be of macromolecule type');
  const uh = new UnitsHandler(column);
  const fnName = uh.getDistanceFunctionName();
  worker.postMessage({values, fnName});
  return new Promise((resolve, reject) => {
    worker.onmessage = ({data: {error, distanceMatrixData}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(distanceMatrixData);
    };
  });
}

// gets index of compressed distance matrix from 2d coordinates
export function dmLinearIndex(size: number) {
  return (i: number, j: number) => {
    const getLinearIndex = (i: number, j: number) => {
      return size * i + j - Math.floor(((i + 2) * (i + 1)) / 2);
    };
    if (i <= j) return getLinearIndex(i, j);
    else return getLinearIndex(j, i);
  };
}
