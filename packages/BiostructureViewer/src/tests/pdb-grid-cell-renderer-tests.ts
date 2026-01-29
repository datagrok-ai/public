import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {before, after, category/*, expect*/, test, testEvent, delay} from '@datagrok-libraries/test/src/test';
import {NglGlServiceBase, getNglGlService} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';

import {awaitGrid} from './utils';
import {IPdbGridCellRenderer} from '../utils/types';
import {PdbGridCellRendererBack} from '../utils/pdb-grid-cell-renderer';

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
    await testEvent(back.onRendered, () => {}, () => {
      grid.invalidate();
    }, 20000 /* prev tests cause timeouts on NGL. "Too many active WebGL contexts" */);
  });
});
