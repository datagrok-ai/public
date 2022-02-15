import * as DG from 'datagrok-api/dg';
import {runKalign} from '../utils/multiple-sequence-alignment';

export async function msaWidget(col: DG.Column): Promise<DG.DataFrame> {
  const msaCol = await runKalign(col, true);
  const table = col.dataFrame;
  table.columns.add(msaCol);
  return table;
}
