// Tests of the multivariate analysis (PLS) results survival in saved projects

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {callMvaTransform, getMvaNames, MvaInput, MvaNames} from '../pls/pls-tools';
import {MVA_MODEL_TAG, MVA_TRANSFORM_FUNC, TITLE} from '../pls/pls-constants';

const DATA_FILE = 'System:AppData/Eda/drugs-props-train.csv';
const COMPONENTS = 3;
const TIMEOUT = 180000;
const CHECK_TIMEOUT = 30000;

const projects: DG.Project[] = [];

/** Open the demo file: the table gets a creation script, which data sync requires */
async function openDataFile(): Promise<DG.TableView> {
  await DG.Func.find({name: 'OpenFile'})[0]
    .prepare({fullPath: DATA_FILE})
    .call(undefined, undefined, {processed: false});

  return grok.shell.tv;
}

/** Run the multivariate analysis the same way the dialog does */
async function runMva(table: DG.DataFrame, componentsOnly: boolean = false): Promise<MvaNames> {
  const numCols = table.columns.toList().filter((col) =>
    ((col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)) &&
    (col.stats.missingValueCount === 0));

  const input: MvaInput = {
    table: table,
    features: numCols.slice(0, -1),
    predict: numCols[numCols.length - 1],
    components: COMPONENTS,
    isQuadratic: false,
    names: undefined,
  };

  const names = getMvaNames(input, componentsOnly);
  await callMvaTransform(input, names, componentsOnly);

  return names;
}

/** Save the table, the analysis tables & the layout to a project, then reopen it */
async function saveAndOpenProject(tv: DG.TableView, names: MvaNames, dataSync: boolean): Promise<void> {
  const project = DG.Project.create();
  project.name = `EDA MVA test (${dataSync ? 'sync' : 'data'})`;
  const tableInfo = tv.dataFrame.getTableInfo();

  if (dataSync) {
    //@ts-ignore: tags is a Dart map proxy
    tableInfo.tags[DG.Tags.DataSync] = 'sync';
    //@ts-ignore
    tableInfo.tags[DG.Tags.CreationScript] = tv.dataFrame.getTag(DG.Tags.CreationScript);
  }

  const layoutInfo = tv.getInfo();
  project.addChild(tableInfo);
  project.addChild(layoutInfo);
  await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
  await grok.dapi.tables.save(tableInfo);

  // the viewers of the analysis tables are docked in the same view, so a project keeps them too
  for (const name of [names.analysisTable, names.explVarTable]) {
    const df = grok.shell.tableByName(name);
    const info = df.getTableInfo();
    project.addChild(info);
    await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.save(info);
  }

  await grok.dapi.views.save(layoutInfo);
  await grok.dapi.projects.save(project);
  projects.push(project);

  grok.shell.closeAll();
  await (await grok.dapi.projects.find(project.id)).open();
}

/** Check that the opened project contains the analysis results */
async function checkResults(tableName: string, names: MvaNames): Promise<void> {
  const expected = names.xScores.concat(names.yScores).concat([names.prediction]);

  await awaitCheck(() => grok.shell.tableNames.includes(tableName),
    `The '${tableName}' table has not been opened`, CHECK_TIMEOUT);

  const table = grok.shell.tableByName(tableName);

  await awaitCheck(() => expected.every((name) => table.col(name) !== null),
    `The analysis columns are missing, expected: ${expected.join(', ')}`, CHECK_TIMEOUT);

  expect(grok.shell.tableNames.includes(names.analysisTable), true,
    `The '${names.analysisTable}' table is missing, opened: ${grok.shell.tableNames.join(', ')}`);
}

category('Multivariate analysis: projects', () => {
  test('Analysis results are added', async () => {
    const tv = await openDataFile();
    const names = await runMva(tv.dataFrame);
    const df = tv.dataFrame;

    for (const name of names.xScores.concat(names.yScores).concat([names.prediction]))
      expect(df.col(name) !== null, true, `The '${name}' column has not been added`);

    const analysis = grok.shell.tableByName(names.analysisTable);

    for (const name of [TITLE.FEATURE, TITLE.REGR_COEFS, `${TITLE.XLOADING}1`, TITLE.VIP])
      expect(analysis.col(name) !== null, true, `The '${name}' column of the analysis table is missing`);

    expect(grok.shell.tableNames.includes(names.explVarTable), true,
      `The '${names.explVarTable}' table has not been created`);

    const model = JSON.parse(df.getTag(MVA_MODEL_TAG)!);
    expect(model.features.length, model.coefficients.length, 'Inconsistent model stored in the table tag');
    expect(model.prediction, names.prediction, 'The model refers to a wrong prediction column');
  }, {timeout: TIMEOUT});

  test('PLS components only', async () => {
    const tv = await openDataFile();
    const names = await runMva(tv.dataFrame, true);

    for (const name of names.xScores)
      expect(tv.dataFrame.col(name) !== null, true, `The '${name}' column has not been added`);

    expect(tv.dataFrame.getTag(MVA_MODEL_TAG) == null, true, 'The model has been stored for PLS components');
  }, {timeout: TIMEOUT});

  test('Analysis is recorded in the creation script', async () => {
    const tv = await openDataFile();
    await runMva(tv.dataFrame);

    const script = tv.dataFrame.getTag(DG.Tags.CreationScript) ?? '';

    expect(script.includes(MVA_TRANSFORM_FUNC), true,
      `The analysis is not in the creation script: ${script}`);
  }, {timeout: TIMEOUT});

  test('Project with data sync restores the results', async () => {
    const tv = await openDataFile();
    const names = await runMva(tv.dataFrame);
    const tableName = tv.dataFrame.name;
    await saveAndOpenProject(tv, names, true);
    await checkResults(tableName, names);
  }, {timeout: TIMEOUT});

  test('Project with uploaded data restores the results', async () => {
    const tv = await openDataFile();
    const names = await runMva(tv.dataFrame);
    const tableName = tv.dataFrame.name;
    await saveAndOpenProject(tv, names, false);
    await checkResults(tableName, names);
  }, {timeout: TIMEOUT});

  after(async () => {
    for (const project of projects)
      await grok.dapi.projects.delete(project);

    projects.length = 0;
    grok.shell.closeAll();
  });
});
