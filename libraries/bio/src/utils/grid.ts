import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export function getGridColByTableCol(grid: DG.Grid, col: DG.Column): DG.GridColumn | null {
  return grid.col(col.name);
}
