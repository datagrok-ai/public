import {category, test, expect, delay, before, testEvent} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import {readDataframe} from './utils';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

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

    async function initSketcher(sw: DG.chem.Sketcher) {
      return new Promise(async (resolve, reject) => {
        sw.sketcherCreated.subscribe(async (_: any) => {
          try {
            resolve(true);
          } catch (error) {
            reject(error);
          }
        });
      });
    }

    async function createFilter(colName: string): Promise<SubstructureFilter> {
      const filter = new SubstructureFilter();
      sketcherDialogs.push(ui.dialog().add(filter.sketcher).show());
      filter.attach(df);
      filter.column = df.col(colName);
      filter.columnName = colName;
      const waitForSketcher = initSketcher(filter.sketcher);
      if (!filter.sketcher.sketcher)
        await waitForSketcher;
      return filter;
    }

    const filter1 = await createFilter('smiles1');
    const filter2 = await createFilter('smiles2');
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
    await testEvent(df.onFilterChanged, (_) => {}, () => { filter1.sketcher.setMolFile(molfile1) }, 7000);
    await testEvent(df.onFilterChanged, (_) => {
      expect(df.filter.trueCount, 2);
      expect(df.filter.get(4), true);
      expect(df.filter.get(5), true);
      expect(df.filter.get(0), false);
    }, () => { filter2.sketcher.setMolFile(molfile2) }, 7000);
    sketcherDialogs.forEach((it) => it.close());
  });
});
  
