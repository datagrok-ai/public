import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Ensures the table view for the given dataframe is active. */
export function checkCurrentView(dataFrame: DG.DataFrame): void {
  if (grok.shell.tv?.dataFrame === dataFrame)
    return;
  const tv = grok.shell.tableView(dataFrame.name);
  if (tv)
    grok.shell.v = tv;
  else
    grok.shell.addTableView(dataFrame);
}
