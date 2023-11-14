import {category, test, expectArray} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {readDataframe} from './utils';
import {scrollTable} from '../utils/demo-utils';
import * as DG from 'datagrok-api/dg';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import { FILTER_SCAFFOLD_TAG, HIGHLIGHT_BY_SCAFFOLD_TAG, SCAFFOLD_TREE_HIGHLIGHT } from '../constants';

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

  test('rdkit grid cell renderer benchmark', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') : await readDataframe('mol1K.csv');
    const rowCount = df.rowCount;
    const col = df.col(DG.Test.isInBenchmark ? 'smiles' : 'molecule')!;
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    const rendndererObj = await renderFunctions[0].apply();
    const moleculeHost = ui.canvas(200, 100);
    
    let start = performance.now();
    //rendering without highlights
    for (let i = 0; i < rowCount; i++) {
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, DG.GridCell.fromValue(col.get(i)));
    }
    console.log(`rendering of ${rowCount} molecules without highlight took ${performance.now() - start} milliseconds`);
   
    col.temp[FILTER_SCAFFOLD_TAG] = 'c1ccccc1';
    col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '[{"color":"#00FF00","molecule":"C1CCCCC1"}]');
    col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '[{"color":"#00FF00","molecule":"C1CCCCC1"}]');
    col.setTag(SCAFFOLD_TREE_HIGHLIGHT, `[{"molecule":"\n  MJ201900                      \n\n  7  7  0  0  0  0  0  0  0  0999 V2000\n    0.1338    0.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5805    0.0107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5805   -0.8143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1338   -1.2268    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8483   -0.8143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8483    0.0107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1338    1.2482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  1  1  0  0  0  0\n  1  7  1  0  0  0  0\nM  END\n","color":"#ffbb78"},{"molecule":"\n  MJ201900                      \n\n  8  8  0  0  0  0  0  0  0  0999 V2000\n    0.1337    0.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5806    0.0107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5806   -0.8144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1337   -1.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8484   -0.8144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8484    0.0107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1337    1.2484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2951    0.4231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  1  1  0  0  0  0\n  1  7  1  0  0  0  0\n  2  8  1  0  0  0  0\nM  END\n","color":"#2ca02c"}]`);

    //rendering with highlights
    start = performance.now();
    for (let i = 0; i < rowCount; i++) {
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, DG.GridCell.fromValue(col.get(i)));
    }
    console.log(`rendering of ${rowCount} molecules with highlight took ${performance.now() - start} milliseconds`);
    
  }, {timeout: 180000});

  test('stereochemistry', async () => {
    const df = await readDataframe('tests/stereochemistry.csv');
    const smiles = df.getCol('smiles').toList();
    const rdKitCellRenderer: RDKitCellRenderer = await grok.functions.call('Chem:rdKitCellRenderer');
    const res = smiles.map((s) => rdKitCellRenderer
      ._fetchMol(s, [], false, false, {}, true).molCtx.useMolBlockWedging);
    expectArray(res, new Array(5).fill(false));
  });
});
