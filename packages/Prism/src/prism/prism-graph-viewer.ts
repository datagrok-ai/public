import * as DG from 'datagrok-api/dg';

import {PrismSheet} from './prism-types';


/** Creates an appropriate Datagrok viewer for a PrismSheet's data. */
export function createGraphViewer(df: DG.DataFrame, sheet: PrismSheet): DG.Viewer | null {
  if (df.rowCount === 0 || df.columns.length < 2)
    return null;

  const format = sheet.dataFormat;
  const hasXColumn = sheet.xTitle != null;
  const hasRowTitles = sheet.rowTitlesPresent;

  // XY data with numeric X → scatter plot
  if (hasXColumn && !hasRowTitles) {
    const xCol = df.columns.byIndex(0);
    if (xCol.type === DG.COLUMN_TYPE.FLOAT || xCol.type === DG.COLUMN_TYPE.INT) {
      const yCol = df.columns.byIndex(1);
      return DG.Viewer.scatterPlot(df, {
        xColumnName: xCol.name,
        yColumnName: yCol.name,
      });
    }
  }

  // Column/grouped data with row titles → bar chart
  if (hasRowTitles || format === 'y_single') {
    const catCol = df.columns.byIndex(0);
    if (catCol.type === DG.COLUMN_TYPE.STRING && df.columns.length >= 2) {
      const valCol = df.columns.byIndex(1);
      return DG.Viewer.barChart(df, {
        splitColumnName: catCol.name,
        valueColumnName: valCol.name,
      });
    }
  }

  // Default: scatter plot with first two numeric columns
  const numCols = df.columns.toList().filter((c) =>
    c.type === DG.COLUMN_TYPE.FLOAT || c.type === DG.COLUMN_TYPE.INT);
  if (numCols.length >= 2) {
    return DG.Viewer.scatterPlot(df, {
      xColumnName: numCols[0].name,
      yColumnName: numCols[1].name,
    });
  }

  return null;
}
