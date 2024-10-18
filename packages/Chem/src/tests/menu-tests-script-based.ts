import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {category, before, expect, test} from '@datagrok-libraries/utils/src/test';
import {molV2000, molV3000, readDataframe} from './utils';


category('top menu script based', () => {
  let spgi100: DG.DataFrame;
  let approvedDrugs100: DG.DataFrame;
  let smiles: DG.DataFrame;

  before(async () => {
    spgi100 = await readDataframe('tests/spgi-100.csv');
    spgi100.rows.removeAt(14, 86);
    approvedDrugs100 = await readDataframe('tests/approved-drugs-100.csv');
    approvedDrugs100.rows.removeAt(14, 86);
    smiles = DG.DataFrame.fromCsv(`Name,smiles
metal_non,CCC(=O)O[Na]
metal_st,CCC(=O)[O-].[Na+]
parent_non,[Na]OC(=O)c1ccccc1
parent_st,O=C([O-])c1ccccc1
norm_non,C[N+](C)=CC=C[O-]
norm_st,CN(C)C=CC=O
reion_non,C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O
reion_st,O=S(O)c1ccc(S(=O)(=O)[O-])cc1
charge_non,O=C([O-])c1ccccc1
charge_st,O=C(O)c1ccccc1
tau_non,C1(=CCCCC1)O
tau_st,O=C1CCCCC1
main_component_non,CCC1=C(C)C=CC(O)=N1.OC(=O)CCC(=O)O
main_component_non_st,CCC1=C(C)C=CC(O)=N1`);
  });

  test('curate.smiles', async () => {
    await curate(DG.Test.isInBenchmark ? grok.data.demo.molecules(1000) : smiles, 'smiles');
  }, {timeout: 60000, benchmark: true});

  test('curate.molV2000', async () => {
    await curate(spgi100, 'Structure');
  }, {timeout: 60000, benchmark: true});

  test('curate.molV3000', async () => {
    await curate(approvedDrugs100, 'molecule');
  }, {timeout: 60000, benchmark: true});

  test('curate.emptyValues', async () => {
    const df = await readDataframe('tests/sar-small_empty_vals.csv');
    await grok.data.detectSemanticTypes(df);
    const t: DG.DataFrame = await grok.functions.call('Chem:Curate', {'data': df, 'molecules': 'smiles',
      'kekulization': true, 'normalization': true, 'reionization': true,
      'neutralization': true, 'tautomerization': true, 'mainFragment': true});
    const col = df.getCol('curated_molecule');
    expect(col.stats.valueCount, 16);
    col.categories.slice(0, -1).forEach((c) => expect(c.includes('C'), true));
  }, {timeout: 60000});

  test('curate.malformedData', async () => {
    const df = await readDataframe('tests/Test_smiles_malformed.csv');
    await grok.data.detectSemanticTypes(df);
    const t: DG.DataFrame = await grok.functions.call('Chem:Curate', {'data': df, 'molecules': 'canonical_smiles',
      'kekulization': true, 'normalization': true, 'reionization': true,
      'neutralization': true, 'tautomerization': true, 'mainFragment': true});
    const col = df.getCol('curated_molecule');
    expect(col.stats.valueCount, 43);
    col.categories.slice(0, -1).forEach((c) => expect(c.includes('C'), true));
  }, {timeout: 60000});

  test('mutate.smiles', async () => {
    await mutate('CN1C(CC(O)C1=O)C1=CN=CC=C1');
  }, {timeout: 60000, benchmark: true});

  test('mutate.molV2000', async () => {
    await mutate(molV2000);
  }, {timeout: 60000, benchmark: true});

  test('mutate.molV3000', async () => {
    await mutate(molV3000);
  }, {timeout: 60000, benchmark: true});

  test('mutate.emptyInput', async () => {
    await mutate('');
  }, {timeout: 60000, benchmark: true});

  test('mutate.malformedInput', async () => {
    await mutate('COc1ccc2c|c(ccc2c1)C(C)C(=O)OCCCc3cccnc3', 0);
  }, {timeout: 60000, benchmark: true});
});


async function curate(df: DG.DataFrame, col: string) {
  const t: DG.DataFrame = await grok.functions.call('Chem:Curate', {'data': df, 'molecules': col,
    'kekulization': true, 'normalization': true, 'reionization': true,
    'neutralization': true, 'tautomerization': true, 'mainFragment': true});
  const cm = df.getCol('curated_molecule');
  if (col !== 'smiles' || DG.Test.isInBenchmark) {
    for (let i = 0; i < t.rowCount; i++) expect(cm.get(i).includes('C'), true);
    return;
  }
  expect(cm.get(0), 'CCC(=O)[O-]');
  expect(cm.get(1), 'CCC(=O)O');
  expect(cm.get(2), 'O=C([O-])C1=CC=CC=C1');
  expect(cm.get(3), 'O=C(O)C1=CC=CC=C1');
  expect(cm.get(4), 'CN(C)C=CC=O');
  expect(cm.get(5), 'CN(C)C=CC=O');
  expect(cm.get(6), 'O=S(O)C1=CC=C(S(=O)(=O)O)C=C1');
  expect(cm.get(7), 'O=S(O)C1=CC=C(S(=O)(=O)O)C=C1');
  expect(cm.get(8), 'O=C(O)C1=CC=CC=C1');
  expect(cm.get(9), 'O=C(O)C1=CC=CC=C1');
  expect(cm.get(10), 'O=C1CCCCC1');
  expect(cm.get(11), 'O=C1CCCCC1');
  expect(cm.get(12), 'CCC1=C(C)C=CC(=O)N1');
  expect(cm.get(13), 'CCC1=C(C)C=CC(=O)N1');
}

async function mutate(molecule: string, expected?: number) {
  const mutations = DG.Test.isInBenchmark && molecule === 'CN1C(CC(O)C1=O)C1=CN=CC=C1' ? 1000 : 10;
  const t: DG.DataFrame = await grok.functions.call('Chem:Mutate', {
    'molecule': molecule,
    'steps': 1,
    'randomize': true,
    'maxRandomResults': mutations,
  });
  expect(t.rowCount, expected ?? mutations);
  const col = t.getCol('mutations');
  for (let i = 0; i < t.rowCount; i++)
    expect(col.get(i).includes('C'), true);
}
