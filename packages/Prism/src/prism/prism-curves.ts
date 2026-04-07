import * as DG from 'datagrok-api/dg';

import {
  IFitChartData,
  IFitChartOptions,
  IFitSeries,
  IFitPoint,
  FIT_FUNCTION_SIGMOID,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';

import {PrismSheet} from './prism-types';


/** Returns true if a sheet contains XY numeric data suitable for curve fitting. */
export function hasXYData(sheet: PrismSheet): boolean {
  if (sheet.xTitle == null && !sheet.rowTitlesPresent)
    return false;
  if (sheet.data.length === 0)
    return false;

  // Check that X values are numeric
  for (const row of sheet.data) {
    if (row.length === 0)
      continue;
    const x = row[0];
    if (x !== '' && isNaN(parseFloat(x)))
      return false;
  }
  return true;
}


/** Converts a PrismSheet with XY data to an IFitChartData object for the Curves plugin. */
export function prismSheetToFitChartData(sheet: PrismSheet): IFitChartData | null {
  if (!hasXYData(sheet))
    return null;

  const xColCount = (sheet.xTitle != null || sheet.rowTitlesPresent) ? 1 : 0;
  if (xColCount === 0)
    return null;

  // Parse X values
  const xValues: (number | null)[] = sheet.data.map((row) => {
    const v = parseFloat(row[0] ?? '');
    return isNaN(v) ? null : v;
  });

  // Determine columns per dataset based on format
  const colsPerDs = getColsPerDataset(sheet);
  const dsCount = sheet.yColumnTitles.length || 1;
  const series: IFitSeries[] = [];

  for (let ds = 0; ds < dsCount; ds++) {
    const startCol = xColCount + ds * colsPerDs;
    const points: IFitPoint[] = [];

    // For replicate formats, collect all replicate values as separate points
    if (sheet.dataFormat === 'y_replicates') {
      const repCount = sheet.replicatesCount;
      for (let ri = 0; ri < sheet.data.length; ri++) {
        const x = xValues[ri];
        if (x == null)
          continue;
        const row = sheet.data[ri];
        for (let rep = 0; rep < repCount; rep++) {
          const ci = startCol + rep;
          if (ci >= row.length)
            continue;
          const y = parseFloat(row[ci]);
          if (!isNaN(y))
            points.push({x, y});
        }
      }
    }
    else {
      // For mean-based formats, use the first (mean) column
      for (let ri = 0; ri < sheet.data.length; ri++) {
        const x = xValues[ri];
        if (x == null)
          continue;
        const row = sheet.data[ri];
        if (startCol >= row.length)
          continue;
        const y = parseFloat(row[startCol]);
        if (!isNaN(y))
          points.push({x, y});
      }
    }

    if (points.length === 0)
      continue;

    series.push({
      name: sheet.yColumnTitles[ds] || `Dataset ${ds + 1}`,
      fitFunction: FIT_FUNCTION_SIGMOID,
      showFitLine: true,
      showPoints: 'points',
      clickToToggle: true,
      droplines: ['IC50'],
      points,
    });
  }

  if (series.length === 0)
    return null;

  const chartOptions: IFitChartOptions = {
    logX: true,
    xAxisName: sheet.xTitle || 'Concentration',
    yAxisName: 'Response',
    title: sheet.title,
  };

  return {chartOptions, series};
}


/** Converts PrismSheets with XY data into a DataFrame with a fit column for Curves rendering. */
export function prismToFitDataFrame(sheets: PrismSheet[]): DG.DataFrame {
  const names: string[] = [];
  const curves: string[] = [];

  for (const sheet of sheets) {
    const chartData = prismSheetToFitChartData(sheet);
    if (!chartData)
      continue;
    names.push(sheet.title || sheet.id);
    curves.push(JSON.stringify(chartData));
  }

  const df = DG.DataFrame.create(names.length);
  df.name = 'Prism Curves';

  const nameCol = df.columns.addNewString('Table');
  for (let i = 0; i < names.length; i++)
    nameCol.set(i, names[i]);

  const curveCol = df.columns.addNewString('Fitted Curve');
  curveCol.semType = 'fit';
  curveCol.meta.cellRenderer = 'fit';
  for (let i = 0; i < curves.length; i++)
    curveCol.set(i, curves[i]);

  return df;
}


/** Returns the number of CSV columns per Y dataset for the given format. */
function getColsPerDataset(sheet: PrismSheet): number {
  switch (sheet.dataFormat) {
    case 'y_single': return 1;
    case 'y_replicates': return sheet.replicatesCount;
    case 'y_sd': case 'y_se': case 'y_cv': return 2;
    case 'y_sd_n': case 'y_se_n': case 'y_cv_n': return 3;
    case 'y_plus_minus': case 'y_high_low': return 3;
    default: return 1;
  }
}
