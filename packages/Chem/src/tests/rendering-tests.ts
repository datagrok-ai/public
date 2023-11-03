import {category, test, expectArray} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {readDataframe} from './utils';
import {scrollTable} from '../utils/demo-utils';
import * as DG from 'datagrok-api/dg';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';

category('rendering', () => {
  test('visual rendering', async () => {
    ui.dialog()
      .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))
      .show();
  });

  test('visual rendering to canvas', async () => {
    const root = ui.div();
    const width = 300;
    const height = 200;
    const canvas1 = ui.canvas(width, height);
    canvas1.id = 'canvas1';
    const canvas2 = ui.canvas(width, height);
    canvas2.id = 'canvas2';
    root.appendChild(canvas1);
    root.appendChild(canvas2);
    // TODO: should be syncronous
    await grok.chem.canvasMol(
      0, 0, width, height, canvas1,
      'COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3',
      'c1ccccc1');
    await grok.chem.canvasMol(
      0, 0, width, height, canvas2,
      'CN1CCC(O)(CC1)c2ccccc2');
    ui
      .dialog({title: 'Test molecule'})
      .add(root)
      .show();
  });

  test('rdkit grid cell renderer', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_1M.zip') :
      await readDataframe('tests/sar-small_test.csv');
    const scrollCycles = DG.Test.isInBenchmark ? 1000 : 10;
    const scrollDelta = DG.Test.isInBenchmark ? 150000 : 10000;
    const tv = grok.shell.addTableView(df);
    const canvas = tv.grid.root.getElementsByTagName('canvas')[2];
    await scrollTable(canvas, scrollDelta, scrollCycles, 5);
  });

  test('stereochemistry', async () => {
    const df = await readDataframe('tests/stereochemistry.csv');
    const smiles = df.getCol('smiles').toList();
    const rdKitCellRenderer: RDKitCellRenderer = await grok.functions.call('Chem:rdKitCellRenderer');
    const res = smiles.map((s) => rdKitCellRenderer
      ._fetchMol(s, [], false, false, {}, true).molCtx.useMolBlockWedging);
    expectArray(res, new Array(5).fill(false));
  });
});
