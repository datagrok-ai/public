import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_oneWayAnovaInWebWorker} from '../../wasm/EDAAPI';
const jerzy = require('jerzy');

// Analysis of Variance (ANOVA)
export async function computeANOVA(table: DG.DataFrame, factor: DG.Column, dependent: DG.Column): Promise<any> {
  const f = factor.toList();
  let d = dependent.toList();
  grok.shell.info('Normal distribution: ' + shapiroWilkTest(d));
  const arrs = transformTable(factor, d);
  const list: DG.Column[] = [];
  // const SWTests = [];
  for (const arr of arrs) {
    const name = arr.shift();
    list.push(DG.Column.fromList(dependent.type, name, arr));
    // SWTests.push(shapiroWilkTest(arr));
  }
  // grok.shell.info('Normal distribution: ' + SWTests.every(Boolean));
  const stdevs = [];
  const df = DG.DataFrame.fromColumns(list);
  for (const col of list)
    stdevs.push(col.stats.stdev);
  grok.shell.info('Equality of variances: ' + (Math.max(...stdevs) < Math.min(...stdevs) * 2));
  if (list.length === 2)
    grok.shell.info('Equality of means: ' + (Math.abs(list[0].stats.avg - list[1].stats.avg) <= 0.05));
  grok.shell.addTableView(df);
  const res: DG.Column = await _oneWayAnovaInWebWorker(df, df.columns);
  const resDf = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'Source', ['Between', 'Within']),
    DG.Column.fromList('double', 'Sum of Squares', [res.get(0), res.get(1)]),
    DG.Column.fromList('int', 'df', [res.get(2), res.get(3)]),
    DG.Column.fromList('double', 'Mean square', [res.get(4), res.get(5)]),
    DG.Column.fromList('double', 'F-test', [res.get(6), null])
  ]);
  grok.shell.addTableView(resDf);
  return 5;
}

function transformTable(factor: DG.Column, values: number[]): any[][] {
  const categories = factor.toList();
  const uniqueCategories = factor.categories;
  const columns: any[] = [];
  const categoryValues: { [key: string]: number[] } = {};
  for (let i = 0; i < uniqueCategories.length; i++) {
    categoryValues[uniqueCategories[i]] = [];
  }
  for (let i = 0; i < categories.length; i++) {
    const category = categories[i];
    const value = values[i];
    if (value)
      categoryValues[category].push(value);
  }
  const min = Math.min(...Object.values(categoryValues).map((v) => v.length));
  Object.keys(categoryValues).forEach((k) => {
    categoryValues[k] = categoryValues[k].slice(0, min);
  });
  for (const category of uniqueCategories) {
    columns.push([category, ...categoryValues[category]]);
  }
  return columns;
}

// Check columns (ANOVA)
export function checkColumns(columns: DG.ColumnList): void {
  if (columns.categorical[Symbol.iterator]().next().value)
    throw new Error('Non numerical columns are selected');

  if (columns.toList().some((col) => col.stats.missingValueCount > 0))
    grok.shell.warning('Columns with missing values are selected, the calculation results may not be accurate.\
      Please exclude columns with missing values from the analysis or fill in the missing data');
}

function shapiroWilkTest(arr: number[], alpha: number = 0.05): boolean {
  const v = new jerzy.Vector(arr.filter(Number));
  return jerzy.Normality.shapiroWilk(v).p >= alpha;
}
