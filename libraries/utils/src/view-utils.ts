import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Ensures the table view for the given dataframe is active. */
export function checkCurrentView(dataFrame: DG.DataFrame): void {
  if (grok.shell.tv?.dataFrame?.dart === dataFrame?.dart)
    return;
  const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame?.dart === dataFrame?.dart);
  if (tv)
    grok.shell.v = tv;
  else
    grok.shell.addTableView(dataFrame);
}
