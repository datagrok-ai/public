import * as DG from 'datagrok-api/dg';
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";

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
