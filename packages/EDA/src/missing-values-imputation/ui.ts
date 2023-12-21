import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { LINK, TITLE, ERROR_MSG } from './constants';

/** Check whether missing values imputation can be applied */
function canImputationBeApplied(col: DG.Column): boolean {
  switch (col.type) {
    case DG.COLUMN_TYPE.INT:
    case DG.COLUMN_TYPE.FLOAT:
    case DG.COLUMN_TYPE.STRING:
    case DG.COLUMN_TYPE.DATE_TIME:
    case DG.COLUMN_TYPE.QNUM:
      return (col.stats.missingValueCount > 0);

    default:
      return false;
  }
}

/** Return array of columns with missing values */
function getColsWithMissingVals(df: DG.DataFrame): DG.Column[] {
  return df.columns.toList().filter((col) => canImputationBeApplied(col));
}

/** */
export function runSimpleImputer(): void {
    
  let df: DG.DataFrame | null = grok.shell.t;
  let cols: DG.Column[];
  let colInput: DG.InputBase<DG.Column> | null = null;
  let strategyInput: DG.InputBase<string> | null = null;
  let inPlaceInput: DG.InputBase<boolean> | null = null;
  let valueInput: DG.InputBase<string | number> | DG.DateInput | null = null;

  const getColInput = (df: DG.DataFrame): DG.InputBase<DG.Column> => {
    ui.columnInput(TITLE.COLUMN, df, null, () => {}, )
  };

  const close

  const dfInput = ui.tableInput(TITLE.TABLE, df, undefined, () => {
    if (dfInput.value === null)
      return;

    cols = getColsWithMissingVals(dfInput.value);

    if (cols.length === 0)
      grok.shell.info(ERROR_MSG.NO_MISSING_VALUES);
    else {
      df = dfInput.value;

      if (colInput === null)
        colInput = getColInput(df);
      else
        colInput.
    }
  });
  
  const dlg = ui.dialog({title: 'Boo', helpUrl: LINK.MISSING_VALUES_IMPUTATION});

  dlg.add(dfInput).onOK(() => {}).show();
}