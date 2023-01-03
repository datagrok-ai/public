import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

/** By Dimitri Petrov */
export function attachDivToGrid(grid: DG.Grid, neighborWidth: number = 100): GridNeighbor {
  // const nRowCount = 100;
  // const nColCount = 5;
  // const dframe: DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);
  // const colDate = dframe.columns.addNewDateTime('New Date');
  // for (let nR = 0; nR < nRowCount; ++nR) {
  //   colDate.set(nR, dayjs());
  // }

  let neighbor: GridNeighbor;
  const eDiv = ui.div();
  // let button = ui.button('Close', () => {
  //   neighbor.close();
  // });
  // eDiv.appendChild(button);
  neighbor = new GridNeighbor(eDiv, grid, neighborWidth);
  return neighbor;
}
