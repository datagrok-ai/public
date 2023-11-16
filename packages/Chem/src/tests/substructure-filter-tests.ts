import {category, test, expect, before, testEvent, delay, expectArray, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import {readDataframe} from './utils';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

const expectedResults: {[key: string]: any} = {
  'oneColumn': [737141248, 593097, 3256025153, 4],
  'malformed': [3029331455, 16127],
  'empty': [524214],
};

category('substructure filters', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('filterBy2Columns', async () => {
    const df = await readDataframe('tests/smiles_2_columns.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles1', df, sketcherDialogs);
    const filter2 = await createFilter('smiles2', df, sketcherDialogs);
    const molfile1 = `
      MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
   -0.1785   -0.1946    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END`;
    const molfile2 = `
      MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
   -0.1785   -0.1946    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
`;

    const terminateFlag1 = 'terminate_substructure_search-tests/smiles_2_columns-smiles1';
    const terminateFlag2 = 'terminate_substructure_search-tests/smiles_2_columns-smiles2';

 
    //finishing first search
    await testEvent(grok.events.onCustomEvent(terminateFlag1), (_) => {},
      () => { filter1.sketcher.setMolFile(molfile1);  }, 7000);

    //finishing fp pre-calculation for filter2
    await testEvent(grok.events.onCustomEvent(terminateFlag2), (_) => {},
      () => {filter2.sketcher.setMolFile(molfile2);}, 7000); 

    //finishing second search
    await testEvent(grok.events.onCustomEvent(terminateFlag2), (_) => {
      expect(df.filter.trueCount, 2);
      expect(df.filter.get(4), true);
      expect(df.filter.get(5), true);
      expect(df.filter.get(0), false);
    }, () => { }, 7000);
    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  });

  test('filterByOneColumn', async () => {
    await testOneColumn('tests/spgi-100.csv', 'Structure', 'c1ccccc1',
      'terminate_substructure_search-tests/spgi-100-Structure', 'oneColumn', 32);
  });

  test('malformed_filterByOneColumn', async () => {
    await testOneColumn('tests/Test_smiles_malformed.csv', 'canonical_smiles', 'c1ccccc1',
      'terminate_substructure_search-tests/Test_smiles_malformed-canonical_smiles', 'malformed', 36);
  });

  test('empty_filterByOneColumn', async () => {
    await testOneColumn('tests/sar-small_empty_vals.csv', 'smiles', 'C1CCCCC1',
      'terminate_substructure_search-tests/sar-small_empty_vals-smiles', 'empty', 16);
  });

  test('terminatedSearch', async () => {
    const df = await readDataframe('tests/smi10K.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter = await createFilter('smiles', df, sketcherDialogs);
    const terminateFlag = 'terminate_substructure_search-tests/smi10K-smiles';
    const substr1 = 'C1CCCCC1';
    const substr2 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter.sketcher.setSmiles(substr1);
    await delay(500);
    //terminating filtering by 1st structure
    await testEvent(grok.events.onCustomEvent(terminateFlag), (_) => {
    }, () => {filter.sketcher.setSmiles(substr2);}, 60000);
    //finishing filtering by 2nd structure
    await testEvent(grok.events.onCustomEvent(terminateFlag), (_) => {
      expect(df.filter.trueCount, 286);
    }, () => {}, 60000);

    sketcherDialogs.forEach((it) => it.close());
    filter.detach();
  }, {timeout: 60000});

  test('filteringMultipleDfs', async () => {
    const df1 = await readDataframe('tests/smi10K.csv');
    const df2 = await readDataframe('tests/smi10K.csv');
    df2.name = 'tests/smi10K (2)';
    await grok.data.detectSemanticTypes(df1);
    await grok.data.detectSemanticTypes(df2);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles', df1, sketcherDialogs);
    const filter2 = await createFilter('smiles', df2, sketcherDialogs);
    const terminateFlag1 = 'terminate_substructure_search-tests/smi10K-smiles';
    const terminateFlag2 = 'terminate_substructure_search-tests/smi10K (2)-smiles';
    const substr1 = 'C1CCCCC1';
    const substr2 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter1.sketcher.setSmiles(substr1);
    await delay(500);
    //opening 2nd filter, sketching a molecule and waiting for 1st filter to complete
    await testEvent(grok.events.onCustomEvent(terminateFlag1), (_) => {
    }, () => {filter2.sketcher.setSmiles(substr2);}, 60000);
    //waiting for 2nd filter to complete
    await testEvent(grok.events.onCustomEvent(terminateFlag2), (_) => {
      expect(df1.filter.trueCount, 462);
      expect(df2.filter.trueCount, 286);
    }, () => {}, 60000);

    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  }, {timeout: 60000});

  test('multipleDfsWithTerminatedSearch', async () => {
    const df1 = await readDataframe('tests/smi10K.csv');
    const df2 = await readDataframe('tests/smi10K.csv');
    df2.name = 'tests/smi10K (2)';
    await grok.data.detectSemanticTypes(df1);
    await grok.data.detectSemanticTypes(df2);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles', df1, sketcherDialogs);
    const filter2 = await createFilter('smiles', df2, sketcherDialogs);
    const terminateFlag1 = 'terminate_substructure_search-tests/smi10K-smiles';
    const terminateFlag2 = 'terminate_substructure_search-tests/smi10K (2)-smiles';
    const substr1 = 'C1CCCCC1';
    const substr2 = 'c1ccccc1';
    const substr3 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter1.sketcher.setSmiles(substr1);
    await delay(500);
    //setting 1st structure to the 2nd filter
    filter2.sketcher.setSmiles(substr2);
    //setting 2nd structure to the 2nd filter and waiting for 1st filter to complete
    await testEvent(grok.events.onCustomEvent(terminateFlag1), (_) => {
    }, () => {filter2.sketcher.setSmiles(substr3);}, 60000);
    //waiting for 2nd filter to complete
    await testEvent(grok.events.onCustomEvent(terminateFlag2), (_) => {
      expect(df1.filter.trueCount, 462);
      expect(df2.filter.trueCount, 286);
    }, () => {}, 60000);

    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  }, {timeout: 60000});
});

async function createFilter(colName: string, df: DG.DataFrame, sketcherDialogs: DG.Dialog[]):
  Promise<SubstructureFilter> {
  const filter = new SubstructureFilter();
  filter.attach(df);
  filter.applyState({columnName: colName});
  sketcherDialogs.push(ui.dialog().add(filter.root).show());
  await ui.tools.waitForElementInDom(filter.sketcher.root);
  filter.column = df.col(colName);
  filter.columnName = colName;
  filter.tableName = df.name;
  await awaitCheck(() => filter.sketcher.sketcher?.isInitialized === true, 'sketcher hasn\'t been initialized', 5000);
  return filter;
}

async function testOneColumn(dfName: string, colName: string, substructure: string, terminateFlag: string,
  expectedKey: string, expectedTrueCount: number) {
  const df = await readDataframe(dfName);
  await grok.data.detectSemanticTypes(df);
  const sketcherDialogs: DG.Dialog[] = [];

  const filter = await createFilter(colName, df, sketcherDialogs);

  filter.sketcher.setSmiles(substructure);
  await testEvent(grok.events.onCustomEvent(terminateFlag), (_) => {
    expect(df.filter.trueCount, expectedTrueCount);
    expectArray(df.filter.getBuffer(), expectedResults[expectedKey]);
  }, () => {}, 7000);
  sketcherDialogs.forEach((it) => it.close());
  filter.detach();
}

