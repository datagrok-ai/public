import * as DG from 'datagrok-api/dg';
import {MwxWorksheet} from './mwx-types';
import {mwxWorksheetToDataFrame} from './mwx-to-dataframe';


/** Creates a TableView for a parsed MWX worksheet. */
export function createMwxTableView(ws: MwxWorksheet): DG.TableView {
  const df = mwxWorksheetToDataFrame(ws);
  return DG.TableView.create(df, false);
}
