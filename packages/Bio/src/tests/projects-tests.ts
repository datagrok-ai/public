import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {readDataframe} from './utils';


category('projects', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  async function createTableView(tableName: string): Promise<DG.TableView> {
    const df = await readDataframe(tableName);
    df.name = tableName.replace('.csv', '');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);
    return view;
  }

  async function saveAndOpenProject(tv: DG.TableView, dataSync?: boolean): Promise<void> {
    const project = DG.Project.create();
    const tableInfo = tv.dataFrame.getTableInfo();
    if (dataSync) {
      //@ts-ignore
      tableInfo.tags[DG.Tags.DataSync] = 'sync';
      //@ts-ignore
      tableInfo.tags[DG.Tags.CreationScript] = grok.shell.tv.dataFrame.getTag(DG.Tags.CreationScript);
    }
    const layoutInfo = tv.getInfo();
    project.addChild(tableInfo);
    project.addChild(layoutInfo);
    await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.views.save(layoutInfo);
    await grok.dapi.projects.save(project);
    const projId = project.id;
    grok.shell.closeAll();
    const p = await grok.dapi.projects.find(projId);
    await p.open();
  }

  async function dataFrameContainsColumns(colArr: string[]): Promise<void> {
    let col = '';
    const getError = () => `${col} hasn't been added to dataframe`;
    await awaitCheck(() => {
      if (!grok.shell.tv.dataFrame)
        return false;
      for (const colName of colArr) {
        if (!grok.shell.tv.dataFrame.col(colName)) {
          col = colName;
          return false;
        }
      }
      return true;
    }, getError(), 5000);
  }

  async function checkViewerAdded(viewerType: string): Promise<void> {
    await awaitCheck(() => {
      for (const v of grok.shell.tv.viewers) {
        if (v.type === viewerType)
          return true;
      }
      return false;
    }, `${viewerType} hasn\'t been added`, 5000);
  }

  async function runSaveAndOpenProjectTest(tableName: string, analysisFunc: (tv: DG.TableView) => Promise<void>,
    colList: string[], viewerType: string, dataSync?: boolean,
    additionalChecks?: (tv: DG.TableView) => Promise<void>) {
    let tv;
    if (dataSync) {
      await DG.Func.find({name: 'OpenFile'})[0].prepare({
        fullPath: `System:AppData/Bio/${tableName}`,
      }).call(undefined, undefined, {processed: false});
      tv = grok.shell.tv;
      await grok.data.detectSemanticTypes(tv.dataFrame);
    } else
      tv = await createTableView(tableName);
    await delay(100);
    await analysisFunc(tv);
    await delay(10);
    await saveAndOpenProject(tv, dataSync);
    await delay(10);
    await dataFrameContainsColumns(colList);
    if (viewerType)
      await checkViewerAdded(viewerType);
    if (additionalChecks)
      await additionalChecks(tv);
  }

  async function runSequenceSpace(tv: DG.TableView): Promise<void> {
    const seqCol = tv.dataFrame.col('sequence')!;
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;
    await DG.Func.find({package: 'Bio', name: 'sequenceSpaceTopMenu'})[0].prepare({
      table: tv.dataFrame,
      molecules: seqCol,
      methodName: 'UMAP',
      similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN,
      plotEmbeddings: true,
      options: {[BYPASS_LARGE_DATA_WARNING]: true},
      clusterEmbeddings: true,
    }).call(undefined, undefined, {processed: false});
    await delay(10);
  }

  async function runActivityCliffs(tv: DG.TableView): Promise<void> {
    const seqCol = tv.dataFrame.col('sequence')!;
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;
    await DG.Func.find({package: 'Bio', name: 'activityCliffs'})[0].prepare({
      table: tv.dataFrame,
      molecules: seqCol,
      activities: tv.dataFrame.col('Activity'),
      similarity: 90,
      methodName: 'UMAP',
      similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN,
      preprocessingFunction: DG.Func.find({name: 'macromoleculePreprocessingFunction', package: 'Bio'})[0],
      options: {[BYPASS_LARGE_DATA_WARNING]: true},
    }).call(undefined, undefined, {processed: false});
    await delay(10);
  }

  async function checkActivityCliffsInit(tv: DG.TableView): Promise<void> {
    let sp: DG.Viewer | null = null;
    for (const v of grok.shell.tv.viewers) {
      if (v.type === DG.VIEWER.SCATTER_PLOT)
        sp = v;
    }
    await awaitCheck(() => {
      const link = sp?.root.getElementsByClassName('scatter_plot_link');
      return !link || !link.length ? false : (link[0] as HTMLElement).innerText.toLowerCase().includes('cliffs');
    }, 'Initialization function hasn\'t been applied on scatter plot', 5000);
  }

  test('sequence_space', async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();

    await runSaveAndOpenProjectTest('tests/100_3_clustests.csv', runSequenceSpace,
      ['sequence', 'Embed_X_1', 'Embed_Y_1', 'Cluster (DBSCAN)'], DG.VIEWER.SCATTER_PLOT);
    await delay(100);

    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  }, {timeout: 60000});

  test('sequence_space_sync', async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();

    await runSaveAndOpenProjectTest('tests/100_3_clustests.csv', runSequenceSpace,
      ['sequence', 'Embed_X_1', 'Embed_Y_1', 'Cluster (DBSCAN)'], DG.VIEWER.SCATTER_PLOT, true);
    await delay(100);

    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  }, {timeout: 60000});

  test('activity_cliffs', async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();

    await runSaveAndOpenProjectTest('tests/100_3_clustests.csv', runActivityCliffs,
      ['sequence', 'Activity', 'Embed_X_1', 'Embed_Y_1'],
      DG.VIEWER.SCATTER_PLOT, false, checkActivityCliffsInit);
    await delay(100);

    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  }, {timeout: 60000});

  test('activity_cliffs_sync', async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();

    await runSaveAndOpenProjectTest('tests/100_3_clustests.csv', runActivityCliffs,
      ['sequence', 'Activity', 'Embed_X_1', 'Embed_Y_1'],
      DG.VIEWER.SCATTER_PLOT, true, checkActivityCliffsInit);
    await delay(100);

    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  }, {timeout: 60000});

  after(async () => {
    grok.shell.closeAll();
  });
});
