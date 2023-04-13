import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import '../css/injected-dendrogram.css';
/** By Dimitri Petrov */
export function attachDivToGrid(grid: DG.Grid, neighborWidth: number = 100): GridNeighbor {
  // const nRowCount = 100;
  // const nColCount = 5;
  // const dframe: DG.DataFrame = grok.data.demo.randomWalk(nRowCount, nColCount);
  // const colDate = dframe.columns.addNewDateTime('New Date');
  // for (let nR = 0; nR < nRowCount; ++nR) {
  //   colDate.set(nR, dayjs());
  // }

  const eDiv = ui.div();
  const button = ui.icons.close(() => {
    neighbor.close();
    //trigger grid rerender
    grid.invalidate();
  }, 'Remove Dendrogram');
  button.classList.add('dendrogram-close-bttn');
  eDiv.appendChild(button);

  const neighbor: GridNeighbor = new GridNeighbor(eDiv, grid, neighborWidth);
  return neighbor;
}
