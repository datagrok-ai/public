import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {before, after, category/*, expect*/, test, testEvent, delay, waitUntil} from '@datagrok-libraries/test/src/test';
import {NglGlServiceBase, getNglGlService} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';

import {awaitGrid} from './utils';
import {IPdbGridCellRenderer} from '../utils/types';
import {PdbGridCellRendererBack, PdbIdGridCellRenderer, pdbIdImageUrl} from '../utils/pdb-grid-cell-renderer';

import {_package} from '../package-test';

category('pdbGridCellRenderer', () => {
  let nglSvc: NglGlServiceBase;

  before(async () => {
    // warm up NglGlDocService
    nglSvc = await getNglGlService();
    await nglSvc.reset();
  });

  after(async () => {
    await nglSvc.reset();
  });

  test('pdbDataCsv', async () => {
    const pdbDataDf: DG.DataFrame =
      await grok.dapi.files.readCsv('System:AppData/BiostructureViewer/pdb_data.csv');
    await pdbDataDf.meta.detectSemanticTypes();
    const pdbDataView: DG.TableView = grok.shell.addTableView(pdbDataDf);

    const grid: DG.Grid = pdbDataView.grid;
    const molGridCol = grid.columns.byName('pdb');
    if (!molGridCol)
      throw new Error(`Column 'pdb' not found.`);

    // const back: IPdbGridCellRenderer = await grok.functions.call(`${_package.name}:getPdbGridCellRenderer`,
    //   {gridCol: molGridCol});
    const gridCell = grid.cell('pdb', 0);
    const back = PdbGridCellRendererBack.getOrCreate(gridCell);

    await delay(0);
    let tim: any = null;
    await testEvent(back.onRendered, () => { tim && clearInterval(tim); }, () => {
      grid.invalidate();
      tim = setInterval(() => grid.invalidate(), 500);
    }, 20000 /* prev tests cause timeouts on NGL. "Too many active WebGL contexts" */);
  });
});

category('pdbIdGridCellRenderer', () => {
  test('image', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('pdb_id', ['1ZP8'])]);
    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const renderer = new PdbIdGridCellRenderer();
    const gridCell = view.grid.cell('pdb_id', 0);
    const canvas = ui.canvas(200, 200);
    renderer.render(canvas.getContext('2d')!, 0, 0, 200, 200, gridCell, gridCell.style);
    await waitUntil(() => renderer.imageCache.has(pdbIdImageUrl('1ZP8')),
      `RCSB has no image of 1ZP8 at ${pdbIdImageUrl('1ZP8')}`, 15000);
  });
});
