import * as DG from 'datagrok-api/dg';
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";


export function getInstalledGridForColumn(colGrid : DG.GridColumn) : DG.Grid | null {
  const dart : any = DG.toDart(colGrid);
  if(!(dart.m_grid instanceof DG.Grid))
    return null;

  return dart.m_grid;
}

export function installGridForColumn(grid : DG.Grid, colGrid : DG.GridColumn) : boolean {
  if(colGrid.grid instanceof DG.Grid)
    return false;

  const dart : any = DG.toDart(colGrid);
  if(dart.m_grid instanceof DG.Grid)
    return false;

  dart.m_grid = grid;
  return true;
}


export function setGridColumnRenderer(colGrid : DG.GridColumn, renderer : GridCellRendererEx) : void {
  const dart : any = DG.toDart(colGrid);
  dart.m_renderer = renderer;
}

export function getGridColumnRenderer(colGrid : DG.GridColumn) : GridCellRendererEx | null {
  const dart : any = DG.toDart(colGrid);
  const renderer = dart.m_renderer;
  if(renderer === undefined)
    return null;

  return renderer;
}
