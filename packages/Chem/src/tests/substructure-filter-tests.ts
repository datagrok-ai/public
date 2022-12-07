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

    const filter1 = new SubstructureFilter();
    const filter2 = new SubstructureFilter();
    filter1.attach(df);
    filter2.attach(df);
    filter1.column = df.col('smiles1');
    filter1.columnName = 'smiles1';
    filter2.column = df.col('smiles2');
    filter2.columnName = 'smiles2';
    
    filter1.sketcher.setMolFile(`
  Ketcher 12 72223422D 1   1.00000     0.00000     0

  5  5  0  0  0  0            999 V2000
    -1.4062   -0.8866    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    -2.1207   -0.4741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.1207    0.3508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6918    0.3508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6918   -0.4741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END`);
    filter2.sketcher.setMolFile(`
  Ketcher 12 72223422D 1   1.00000     0.00000     0

  5  5  0  0  0  0            999 V2000
    -1.4062   -0.8866    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    -2.1207   -0.4741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.1207    0.3508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6918    0.3508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6918   -0.4741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END
`); 
    await testEvent(df.onFilterChanged, (_) => {}, () => { filter1._onSketchChanged() }, 5000);
    await testEvent(df.onFilterChanged, (_) => {
      expect(df.filter.trueCount, 2);
      expect(df.filter.get(4), true);
      expect(df.filter.get(5), true);
      expect(df.filter.get(0), false);
    }, () => { filter2._onSketchChanged() }, 5000);
  });
});
  
