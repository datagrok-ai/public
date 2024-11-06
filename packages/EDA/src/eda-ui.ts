// Custom UI for Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Rename PCA columns
export function renamePCAcolumns(pcaTable: DG.DataFrame): DG.DataFrame {
  for (const col of pcaTable.columns.toList())
    col.name = 'PCA' + col.name;

  return pcaTable;
}

// Adds prefix to each column name
export function addPrefixToEachColumnName(prefix: string, columns: DG.ColumnList): void {
  for (const col of columns.toList())
    col.name = prefix + col.name;
}
