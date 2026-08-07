import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {changeCurvesOptions} from '../fit/fit-options';
import {getColumnChartOptions, getOrCreateParsedChartData} from '../fit/fit-chart-data';

const DEMO_FILE = 'curves-demo.csv';
// every cell of this column declares showCurveConfidenceInterval, so a column-level "off" only shows
// up if the option outranks the data rather than being written into it
const CURVE_COLUMN = 'multiple styled series';

async function createTableView(tableName: string): Promise<DG.TableView> {
  const df = await grok.dapi.files.readCsv(`System:AppData/Curves/${tableName}`);
  df.name = tableName.replace('.csv', '');
  await grok.data.detectSemanticTypes(df);
  return grok.shell.addTableView(df);
}

async function saveAndOpenProject(tv: DG.TableView): Promise<void> {
  const project = DG.Project.create();
  project.name = 'Curves test project';
  const tableInfo = tv.dataFrame.getTableInfo();
  const layoutInfo = tv.getInfo();
  project.addChild(tableInfo);
  project.addChild(layoutInfo);
  await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.views.save(layoutInfo);
  await grok.dapi.projects.save(project);
  const projId = project.id;
  // the fit cell renderer keeps painting while the view tears down, and closeAll() surfaces that as
  // an unhandled exception - let the grid settle first
  await delay(500);
  grok.shell.closeAll();
  await delay(500);
  const p = await grok.dapi.projects.find(projId);
  await p.open();
}

async function dataFrameContainsColumn(colName: string): Promise<void> {
  await awaitCheck(() => grok.shell.tv?.dataFrame?.col(colName) != null,
    `${colName} is missing from the reopened project`, 30000);
}

function turnOffConfidenceIntervals(tv: DG.TableView): void {
  if (tv.dataFrame.col(CURVE_COLUMN) === null)
    throw new Error(`no ${CURVE_COLUMN} column - the table has: ${tv.dataFrame.columns.names().join(', ')}`);
  changeCurvesOptions(tv.grid.cell(CURVE_COLUMN, 0),
    {property: {name: 'showCurveConfidenceInterval'}, value: false} as unknown as DG.InputBase,
    'seriesOptions', 'Column');
}

async function runSaveAndOpenProjectTest(): Promise<void> {
  const tv = await createTableView(DEMO_FILE);
  await delay(100);
  turnOffConfidenceIntervals(tv);
  await delay(10);
  await saveAndOpenProject(tv);
  await delay(10);
  await dataFrameContainsColumn(CURVE_COLUMN);
  // the tags are restored after the table, so wait for the option rather than for the column
  await awaitCheck(() => getColumnChartOptions(grok.shell.tv.dataFrame.col(CURVE_COLUMN)!)
    .seriesOptions?.showCurveConfidenceInterval === false,
  'the column tag did not come back with the project', 30000);

  const df = grok.shell.tv.dataFrame;
  expect(getOrCreateParsedChartData(df.cell(0, CURVE_COLUMN)).series![0].showCurveConfidenceInterval, false,
    'the column-level option did not outrank the data after reopening');
}

category('projects', () => {
  test('column-level option survives a project', async () => {
    await runSaveAndOpenProjectTest();
    await delay(100);
  }, {timeout: 120000});
});
