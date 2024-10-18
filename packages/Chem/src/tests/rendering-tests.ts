import {category, test, expectArray, before} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {readDataframe} from './utils';
import {scrollTable} from '../utils/demo-utils';
import * as DG from 'datagrok-api/dg';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import { FILTER_SCAFFOLD_TAG, HIGHLIGHT_BY_SCAFFOLD_TAG, SCAFFOLD_TREE_HIGHLIGHT } from '../constants';
import { convertMolNotation } from '../package';
import { _convertMolNotation } from '../utils/convert-notation-utils';

category('rendering', () => {
  
  let rdkitModule: any;
  
  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
  });

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

  test('rdkit grid cell renderer scroll', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('System:DemoFiles/chem/smiles_1M.zip') :
      await readDataframe('tests/sar-small_test.csv');
    const scrollCycles = DG.Test.isInBenchmark ? 1000 : 10;
    const scrollDelta = DG.Test.isInBenchmark ? 150000 : 10000;
    const tv = grok.shell.addTableView(df);
    const canvas = tv.grid.root.getElementsByTagName('canvas')[2];
    await scrollTable(canvas, scrollDelta, scrollCycles, 5);
  }, {benchmark: true});

  test('rdkit grid cell renderer', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') : await readDataframe('mol1K.csv');
    const grid = df.plot.grid();
    const rowCount = df.rowCount;
    const colName = DG.Test.isInBenchmark ? 'smiles' : 'molecule'
    const col = df.col(colName)!;
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    const rendndererObj = await renderFunctions[0].apply();
    const moleculeHost = ui.canvas(200, 100);
    
    let start = performance.now();
    //rendering without highlights
    for (let i = 0; i < rowCount; i++) {
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, grid.cell(colName, i));
    }
    console.log(`rendering of ${rowCount} molecules without highlight took ${performance.now() - start} milliseconds`);
    
  }, {timeout: 180000, benchmark: true});

  test('rdkit grid cell renderer with highlights', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') : await readDataframe('mol1K.csv');
    const grid = df.plot.grid();
    const rowCount = df.rowCount;
    const colName = DG.Test.isInBenchmark ? 'smiles' : 'molecule'
    const col = df.col(colName)!;
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    const rendndererObj = await renderFunctions[0].apply();
    const moleculeHost = ui.canvas(200, 100);
    
    let start = performance.now();
   
    col.temp[FILTER_SCAFFOLD_TAG] = JSON.stringify([{
      molecule: _convertMolNotation('c1ccccc1', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
      isSuperstructure: false
    }]);
    col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '[{"color":"#00FF00","molecule":"C1CCCCC1"}]');
    const scaffoldTag = JSON.stringify([
      {
        molecule: _convertMolNotation('Cc1ccccc1', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
        color: '#ffbb78'
      },
      {
        molecule: _convertMolNotation('Cc1ccccc1C', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
        color: '#2ca02c'
      },
    ]);
    col.setTag(SCAFFOLD_TREE_HIGHLIGHT, scaffoldTag);
    //rendering with highlights
    start = performance.now();
    for (let i = 0; i < rowCount; i++) {
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, grid.cell(colName, i));
    }
    console.log(`rendering of ${rowCount} molecules with highlight took ${performance.now() - start} milliseconds`);
    
  }, {timeout: 180000, benchmark: true});

  test('stereochemistry', async () => {
    const df = await readDataframe('tests/stereochemistry.csv');
    const smiles = df.getCol('smiles').toList();
    const rdKitCellRenderer: RDKitCellRenderer = await grok.functions.call('Chem:rdKitCellRenderer');
    const res = smiles.map((s) => rdKitCellRenderer
      ._fetchMol(s, [], false, false, {}, true).molCtx.useMolBlockWedging);
    expectArray(res, new Array(5).fill(false));
  });
});
