import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {category, test, before, after, expect, expectArray} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {runStructuralAlertsDetection} from '../panels/structural-alerts';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {elementalAnalysis} from '../../src/package';
import {readDataframe} from './utils';


category('screening tools', () => {
  let spgi100: DG.DataFrame;
  let approvedDrugs100: DG.DataFrame;
  let molecules: DG.DataFrame;

  before(async () => {
    chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
    chemCommonRdKit.initRdKitModuleLocal();
    spgi100 = await readDataframe('tests/spgi-100.csv');
    approvedDrugs100 = await readDataframe('tests/approved-drugs-100.csv');
    molecules = grok.data.demo.molecules(100);
    await grok.data.detectSemanticTypes(spgi100);
    await grok.data.detectSemanticTypes(approvedDrugs100);
    await grok.data.detectSemanticTypes(molecules);
  });

  test('elementalAnalysis.smiles', async () => {
    let df: DG.DataFrame;
    if (DG.Test.isInBenchmark)
      df = await readDataframe('smiles.csv');
    else 
      df = molecules.clone();
    const tv = grok.shell.addTableView(df);
    elementalAnalysis(df, df.getCol(DG.Test.isInBenchmark ? 'canonical_smiles' : 'smiles'), false, false);
    tv.close();
    expect(df.columns.length, DG.Test.isInBenchmark ? 20 : 11);
  });

  test('elementalAnalysis.molV2000', async () => {
    const df = spgi100.clone();
    elementalAnalysis(df, df.getCol('Structure'), false, false);
    expect(df.columns.length, 95);
  });

  test('elementalAnalysis.molV3000', async () => {
    const df = approvedDrugs100.clone();
    elementalAnalysis(df, df.getCol('molecule'), false, false);
    expect(df.columns.length, 41);
  });

  test('elementalAnalysis.emptyValues', async () => {
    const df = await readDataframe('tests/sar-small_empty_vals.csv');
    await grok.data.detectSemanticTypes(df);
    elementalAnalysis(df, df.getCol('smiles'), false, false);
    expect(df.columns.length, 6);
    expectArray(Array.from(df.row(0).cells).map((c) => c.value), ['', 0, 0, 0, 0, 0]);
  }, {skipReason: 'GROK-12227'});

  test('elementalAnalysis.malformedData', async () => {
    const df = await readDataframe('tests/Test_smiles_malformed.csv');
    await grok.data.detectSemanticTypes(df);
    elementalAnalysis(df, df.getCol('canonical_smiles'), false, false);
    expect(df.columns.length, 29);
    expect(Array.from(df.row(40).cells).map((c) => c.value).join(''),
      '1480010COc1ccc2cc(ccc2c1)C(C)C(=O)OC|CCc3cccnc300040203710400.272729992866516126340');
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
});


// To do: move to separate Benchmarks category
category('screening tools: benchmarks', () => {
  before(async () => {
    chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
    chemCommonRdKit.initRdKitModuleLocal();
  });

  test('structural alerts', async () => {
    const alertsDf = DG.DataFrame.fromCsv(await _package.files.readAsText('alert-collection.csv'));
    const ruleSetCol = alertsDf.getCol('rule_set_name');
    const smartsCol = alertsDf.getCol('smarts');
    const ruleIdCol = alertsDf.getCol('rule_id');
    const rdkitModule = chemCommonRdKit.getRdKitModule();

    const smartsMap = new Map<string, RDMol>();
    for (let i = 0; i < alertsDf.rowCount; i++)
      smartsMap.set(ruleIdCol.get(i), rdkitModule.get_qmol(smartsCol.get(i)));

    const sarSmall = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/smi10K.csv'));
    const smilesCol = sarSmall.getCol('smiles');
    const ruleSetList = ['BMS', 'Dundee', 'Glaxo', 'Inpharmatica', 'LINT', 'MLSMR', 'PAINS', 'SureChEMBL'];

    DG.time('Structural Alerts', () => {
      runStructuralAlertsDetection(sarSmall, ruleSetList, smilesCol, ruleSetCol, ruleIdCol, smartsMap, rdkitModule);
    });
  }, {skipReason: '#1193'});

  test('elementalAnalysis', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(await _package.files.readAsText('test.csv'));
    await grok.data.detectSemanticTypes(df);
    const col: DG.Column = df.getCol('molecule');
    DG.time('Elemental Analysis', async () => {
      await elementalAnalysis(df, col, false, false);
    });
  });
});
