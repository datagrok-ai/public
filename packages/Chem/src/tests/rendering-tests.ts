import {category, test, expect, expectArray, before, awaitCheck} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {readDataframe} from './utils';
import {scrollTable} from '../utils/demo-utils';
import * as DG from 'datagrok-api/dg';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import {showEnlargedMolecule} from '../rendering/render-molecule';
import {FILTER_SCAFFOLD_TAG, HIGHLIGHT_BY_SCAFFOLD_TAG, SCAFFOLD_TREE_HIGHLIGHT} from '../constants';
import {_convertMolNotation} from '../utils/convert-notation-utils';

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
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
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
    const colName = DG.Test.isInBenchmark ? 'smiles' : 'molecule';
    const col = df.col(colName)!;
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    const rendndererObj = await renderFunctions[0].apply();
    const moleculeHost = ui.canvas(200, 100);

    const start = performance.now();
    //rendering without highlights
    for (let i = 0; i < rowCount; i++)
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, grid.cell(colName, i));

    console.log(`rendering of ${rowCount} molecules without highlight took ${performance.now() - start} milliseconds`);
  }, {timeout: 180000, benchmark: true});

  test('rdkit grid cell renderer with highlights', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') : await readDataframe('mol1K.csv');
    const grid = df.plot.grid();
    const rowCount = df.rowCount;
    const colName = DG.Test.isInBenchmark ? 'smiles' : 'molecule';
    const col = df.col(colName)!;
    const renderFunctions = DG.Func.find({meta: {chemRendererName: 'RDKit'}});
    const rendndererObj = await renderFunctions[0].apply();
    const moleculeHost = ui.canvas(200, 100);

    let start = performance.now();

    col.temp[FILTER_SCAFFOLD_TAG] = JSON.stringify([{
      molecule: _convertMolNotation('c1ccccc1', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
      isSuperstructure: false,
    }]);
    col.setTag(HIGHLIGHT_BY_SCAFFOLD_TAG, '[{"color":"#00FF00","molecule":"C1CCCCC1"}]');
    const scaffoldTag = JSON.stringify([
      {
        molecule: _convertMolNotation('Cc1ccccc1', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
        color: '#ffbb78',
      },
      {
        molecule: _convertMolNotation('Cc1ccccc1C', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule),
        color: '#2ca02c',
      },
    ]);
    col.setTag(SCAFFOLD_TREE_HIGHLIGHT, scaffoldTag);
    //rendering with highlights
    start = performance.now();
    for (let i = 0; i < rowCount; i++)
      rendndererObj.render(moleculeHost.getContext('2d')!, 0, 0, 200, 100, grid.cell(colName, i));

    console.log(`rendering of ${rowCount} molecules with highlight took ${performance.now() - start} milliseconds`);
  }, {timeout: 180000, benchmark: true});

  test('mixture cell renderer', async () => {
    const df = await readDataframe('test_mixtures.csv');
    await grok.data.detectSemanticTypes(df);
    expect(df.col('mixture')!.semType, 'ChemicalMixture');
    const grid = df.plot.grid();
    const renderer: DG.GridCellRenderer =
      await DG.Func.find({package: 'Chem', meta: {cellType: 'ChemicalMixture', role: 'cellRenderer'}})[0].apply();
    const r = window.devicePixelRatio;
    // the renderer works in device pixels: a fixed device-pixel canvas gets the same picture at any devicePixelRatio
    const [w, h, margin] = [600, 200, 4];
    const leafComponents = [2, 3, 3, 3];
    for (let row = 0; row < df.rowCount; row++) {
      const canvas = ui.canvas(w, h);
      const g = canvas.getContext('2d')!;
      const cell = grid.cell('mixture', row);
      renderer.render(g, 0, 0, w / r, h / r, cell, cell.style);
      const pixels = g.getImageData(0, 0, w, h).data;
      const slot = (w - 2 * margin) / leafComponents[row];
      const ink = new Array(leafComponents[row]).fill(0);
      for (let y = 0; y < h; y++) {
        for (let x = margin; x < w - margin; x++) {
          const i = (y * w + x) * 4;
          if (pixels[i + 3] < 128 || Math.min(pixels[i], pixels[i + 1], pixels[i + 2]) > 230)
            continue;
          ink[Math.min(Math.floor((x - margin) / slot), ink.length - 1)]++;
        }
      }
      expect(ink.every((n) => n > 10), true, `row ${row}: a component drew no ink (${ink.join(', ')})`);
    }
  });

  test('stereochemistry', async () => {
    const df = await readDataframe('tests/stereochemistry.csv');
    const smiles = df.getCol('smiles').toList();
    const rdKitCellRenderer: RDKitCellRenderer = await grok.functions.call('Chem:rdKitCellRenderer');
    const res = smiles.map((s) => rdKitCellRenderer
      ._fetchMol(s, [], false, false, {}, true).molCtx.useMolBlockWedging);
    expectArray(res, new Array(5).fill(false));
  });

  test('Enlarge draws the molecule again at the window\'s size', async () => {
    const dialog = showEnlargedMolecule('O=C(O)c1ccccc1');
    try {
      const canvas = dialog.root.querySelector('canvas')!;
      expect(canvas.style.width, `${Math.max(200, window.innerWidth - 200)}px`);
      await awaitCheck(() => {
        const data = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height).data;
        return data.some((v, i) => i % 4 === 0 && v < 100 && data[i + 3] > 200); // a bond: dark and opaque
      }, 'nothing drawn', 5000);
    } finally {
      dialog.close();
    }
  });
});
