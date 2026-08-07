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

async function saveAndOpenProject(tv: DG.TableView, dataSync?: boolean): Promise<void> {
  const project = DG.Project.create();
  project.name = 'Curves test project';
  const tableInfo = tv.dataFrame.getTableInfo();
  if (dataSync) {
    //@ts-ignore
    tableInfo.tags[DG.Tags.DataSync] = 'sync';
    //@ts-ignore
    tableInfo.tags[DG.Tags.CreationScript] = tv.dataFrame.getTag(DG.Tags.CreationScript);
  }
  // saveLayout, not getInfo: only saveLayout runs saveSyncTags, which is what puts the column tags
  // into the layout - without them a refetched table comes back with nothing but its own data
  const layoutInfo = tv.saveLayout();
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

async function runSaveAndOpenProjectTest(dataSync?: boolean): Promise<void> {
  let tv: DG.TableView;
  if (dataSync) {
    await DG.Func.find({name: 'OpenFile'})[0].prepare({
      fullPath: `System:AppData/Curves/${DEMO_FILE}`,
    }).call(undefined, undefined, {processed: false});
    tv = grok.shell.tv;
    await grok.data.detectSemanticTypes(tv.dataFrame);
  } else { tv = await createTableView(DEMO_FILE); }
  await delay(100);
  turnOffConfidenceIntervals(tv);
  await delay(10);
  await saveAndOpenProject(tv, dataSync);
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

  test('column-level option survives a datasync project', async () => {
    // the table is rebuilt from its creation script, so the cells come back declaring the interval
    // and only the tag the layout carries can turn it off
    await runSaveAndOpenProjectTest(true);
    await delay(100);
  }, {timeout: 120000});
});
