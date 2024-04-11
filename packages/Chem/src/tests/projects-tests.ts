import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {after, awaitCheck, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import { createTableView } from './utils';


category('projects', () => {
  test('r-group-analysis', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runRGroupAnalysis,
      ['smiles', 'Core', 'R1', 'R2', 'R3', 'isHit'], DG.VIEWER.TRELLIS_PLOT);
  });

  test('inchi', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchi,
      ['smiles', 'inchi'], '');
  });

  test('inchi_key', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchiKeys,
      ['smiles', 'inchi_key'], '');
  });

  after(async () => {
    grok.shell.closeAll();
  });
});

//colList: list of column which have been added to dataframe by analysis
//viewerType: viewer which has been created by analysis
async function runSaveAndOpenProjectTest(tableName: string, analysisFunc: (tv: DG.TableView) => Promise<void>,
colList: string[], viewerType: string, dataSync?: boolean) {
  const tv = await createTableView(tableName);
  await analysisFunc(tv);
  await delay(10);
  await saveAndOpenProject(tv);
  await dataFrameContainsColumns(colList);
  if (viewerType)
    await checkViewerAdded(viewerType);
}

async function runToInchiKeys(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:addInchisKeysTopMenu`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
  });
}

async function runToInchi(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:addInchisTopMenu`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
  });
}

async function runRGroupAnalysis(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:rGroupDecomposition`, {
    df: tv.dataFrame,
    molColName: 'smiles',
    core: '[#8]=[#6]1-[#6]-[#7]=[#6](-[#6])-[#6]2:[#6](-[#7]-1):[#6]:[#6]:[#6]:[#6]:2',
    rGroupName: 'R',
    rGroupChunkSize: '5',
    rGroupMatchingStrategy: 'Greedy',
    rGroupAlignment: 'MCS',
    visualAnalysis: true
  });
  tv.trellisPlot({
    xColumnNames: ['R1'],
    yColumnNames: ['R2'],
  });
}

async function saveAndOpenProject(tv: DG.TableView): Promise<void> {
  let project = DG.Project.create();
  let tableInfo = tv.dataFrame.getTableInfo();
  let layoutInfo = tv.getInfo();
  project.addChild(tableInfo);
  project.addChild(layoutInfo);
  await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.views.save(layoutInfo);
  await grok.dapi.projects.save(project);
  const projId = project.id;
  grok.shell.closeAll();
  const p = await grok.dapi.projects.find(projId);
  p.open();
}

async function dataFrameContainsColumns(colArr: string[]): Promise<void> {
  let col = '';
  await awaitCheck(() => {
    if (!grok.shell.tv.dataFrame)
      return false;
    for (const colName of colArr) {
      if (!grok.shell.tv.dataFrame.columns.names().includes(colName)) {
        col = colName;
        return false;
      }
    }
    return true;
  }, `${col} hasn't been added to dataframe`, 5000);
}

async function checkViewerAdded(viewerType: string): Promise<void> {
  await awaitCheck(() => {
    for (let v of grok.shell.tv.viewers) {
      if (v.type === viewerType)
        return true;
    }
    return false;
  }, `${viewerType} hasn\'t been added`, 5000);
}

