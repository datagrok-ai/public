import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {after, awaitCheck, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import { createTableView } from './utils';
import { RGroupDecompRes } from '../analysis/r-group-analysis';
import { _package } from '../package-test';


category('projects', () => {
  test('r-group-analysis', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runRGroupAnalysis,
      ['smiles', 'Core', 'R1', 'R2', 'R3', 'isHit'], DG.VIEWER.TRELLIS_PLOT);
  });

  // test('r-group-analysis-sync', async () => {
  //   await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runRGroupAnalysis,
  //     ['smiles', 'Core', 'R1', 'R2', 'R3', 'isHit'], DG.VIEWER.TRELLIS_PLOT, true);
  // });

  test('inchi', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchi, ['smiles', 'inchi'], '');
  });

  test('inchi_sync', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchi, ['smiles', 'inchi'], '', true);
  });

  test('inchi_key', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchiKeys, ['smiles', 'inchi_key'], '');
  });

  // test('inchi_key_sync', async () => {
  //   await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runToInchiKeys, ['smiles', 'inchi_key'], '', true);
  // });

  test('toxicity_risks', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runAddChemRisksColumns,
      ['smiles', 'Mutagenicity'], '');
  });

  test('properties', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runAddChemPropertiesColumns,
      ['smiles', 'MW'], '');
  });

  test('curate', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runCurate,
      ['smiles', 'curated_molecule'], '');
  });

  test('names_to_smiles', async () => {
    await runSaveAndOpenProjectTest('tests/names_to_smiles.csv', runNamesToSmiles,
      ['Name', 'canonical_smiles'], '');
  });

  test('convert_notation', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runConvertNotation,
      ['smiles', 'smiles_molblock'], '');
  });

  test('elemental_analysis', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runElementalAnalysis,
      ['smiles', 'C', 'N', 'O', 'Cl', 'Molecule Charge'], DG.VIEWER.RADAR_VIEWER);
    expect(grok.shell.tv.grid.col('elements (smiles)') != null);
  });

  test('structural_alerts', async () => {
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runStructuralAlerts,
    ['smiles', 'PAINS (smiles)'], '');
  });

  test('chemical_space', async () => {
    //column '~smiles.Morgan' is not saved to project since it is an object type column
    await runSaveAndOpenProjectTest('tests/sar-small_test.csv', runChemicalSpace,
    ['smiles', 'Embed_X_1', 'Embed_Y_1', 'Cluster (DBSCAN)'], DG.VIEWER.SCATTER_PLOT);
    //need delay to avoid unhandled exception when calling closeAll()
    await delay(100);
  });

  after(async () => {
    grok.shell.closeAll();
  });
});

//colList: list of column which have been added to dataframe by analysis
//viewerType: viewer which has been created by analysis
async function runSaveAndOpenProjectTest(tableName: string, analysisFunc: (tv: DG.TableView) => Promise<void>,
colList: string[], viewerType: string, dataSync?: boolean) {
  //const tv = await createTableView(tableName);

  const df = await DG.Func.find({ name: 'OpenFile'})[0].prepare({
    fullPath: `System:AppData/Chem/${tableName}`
  }).call(undefined, undefined, { processed: false });

  const tv = grok.shell.tv;
  await delay(500);
  await analysisFunc(tv);
  await delay(10);
  await saveAndOpenProject(tv, dataSync);
  await dataFrameContainsColumns(colList);
  if (viewerType)
    await checkViewerAdded(viewerType);
}

async function runChemicalSpace(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:chemSpaceTopMenu`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    methodName: 'UMAP',
    similarityMetric: 'Tanimoto',
    plotEmbeddings: true,
    options: undefined,
    preprocessingFunction: undefined,
    clusterEmbeddings: true
  });
  //need for scatter plot to render
  await delay(10);
}

async function runStructuralAlerts(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:runStructuralAlerts`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    pains: true,
    bms: false,
    sureChembl: false,
    mlsmr: false,
    dandee: false,
    inpharmatica: false,
    lint: false,
    glaxo: false
  });
}

async function runElementalAnalysis(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:elementalAnalysis`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    radarViewer: true,
    radarGrid: true
  });
}

async function runConvertNotation(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:convertNotation`, {
    data: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    targetNotation: 'molblock',
    overwrite: false
  });
}

async function runNamesToSmiles(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:namesToSmiles`, {
    data: tv.dataFrame,
    names: tv.dataFrame.col('Name'),
  });
}

async function runCurate(tv: DG.TableView): Promise<void> {
  const df = await grok.functions.call(`Chem:Curate`, {
    data: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    kekulization: true,
    normalization: false,
    reionization: false,
    neutralization: false,
    tautomerization: false,
    mainFragment: false,
  });
  tv.dataFrame.columns.add(df.col('curated_molecule'));
}

async function runAddChemPropertiesColumns(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:addChemPropertiesColumns`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    MW: true,
    HBA: false,
    HBD: false,
    logP: false,
    logS: false,
    PSA: false,
    rotatableBonds: false,
    stereoCenters: false,
    moleculeCharge: false,
  });
}

async function runToInchiKeys(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:addInchisKeysTopMenu`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
  });
}

async function runToInchi(tv: DG.TableView): Promise<void> {
  await DG.Func.find({ package: 'Chem', name: 'addInchisTopMenu' })[0].prepare({
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
  }).call(undefined, undefined, { processed: false });
}

async function runAddChemRisksColumns(tv: DG.TableView): Promise<void> {
  await grok.functions.call(`Chem:addChemRisksColumns`, {
    table: tv.dataFrame,
    molecules: tv.dataFrame.col('smiles'),
    mutagenicity: true,
    tumorigenicity: false,
    irritatingEffects: false,
    reproductiveEffects: false,
  });
}

async function runRGroupAnalysis(tv: DG.TableView): Promise<void> {
  const funcCall = await DG.Func.find({ package: 'Chem', name: 'rGroupDecomposition' })[0].prepare({
    df: tv.dataFrame,
    molColName: 'smiles',
    core: '[#8]=[#6]1-[#6]-[#7]=[#6](-[#6])-[#6]2:[#6](-[#7]-1):[#6]:[#6]:[#6]:[#6]:2',
    rGroupName: 'R',
    rGroupChunkSize: '5',
    rGroupMatchingStrategy: 'Greedy',
    rGroupAlignment: 'MCS',
    visualAnalysis: true
  }).call(undefined, undefined, { processed: false });
  const res: RGroupDecompRes = funcCall.getOutputParamValue();
  tv.trellisPlot({
    xColumnNames: res.xAxisColName,
    yColumnNames: res.yAxisColName,
  });
}

async function saveAndOpenProject(tv: DG.TableView, dataSync?: boolean): Promise<void> {
  let project = DG.Project.create();
  let tableInfo = tv.dataFrame.getTableInfo();
  if (dataSync)
    tv.dataFrame.setTag('.data-sync', 'sync');
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
    for (let v of grok.shell.tv.viewers) {
      if (v.type === viewerType)
        return true;
    }
    return false;
  }, `${viewerType} hasn\'t been added`, 5000);
}

