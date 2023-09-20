import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {DistanceMatrixService} from '@datagrok-libraries/ml/src/distance-matrix';

export async function calcMmDistanceMatrix(column: DG.Column<any>): Promise<Float32Array> {
  const values = column.toList();
  if (column.semType !== DG.SEMTYPE.MACROMOLECULE)
    throw new Error('Column has to be of macromolecule type');
  const uh = UnitsHandler.getOrCreate(column);
  const fnName = uh.getDistanceFunctionName();
  const distanceMatrixService = new DistanceMatrixService(true, false);
  const dm = await distanceMatrixService.calc(values, fnName);
  distanceMatrixService.terminate();
  return dm;
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
