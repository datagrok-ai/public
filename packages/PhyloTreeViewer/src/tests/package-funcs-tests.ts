import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {
  after,
  before,
  category,
  test,
  expect,
  expectArray,
  expectObject,
  waitUntil
} from '@datagrok-libraries/test/src/test';

import {PhylocanvasGlService} from '../utils/phylocanvas-gl-service';
import {TreeToGridApp} from '../apps/tree-to-grid-app';
import {awaitGrid} from './utils';

import {_package} from '../package-test';

/** Tests for package functions, test apps, file previews, file handlers, ... */
category('packageFuncs', () => {
  test('phylocanvasGlViewerApp', async () => {
    await grok.functions.call(`${_package.name}:phylocanvasGlViewerApp`, {});
    await waitUntil(() => true, 'Error');
  });

  test('injectTreeToGrid', async () => {
    const nwk = await grok.dapi.files.readAsText(`System:AppData/${_package.name}/data/tree95.nwk`);
    const df: DG.DataFrame = await grok.dapi.files.readCsv(`System.AppData/${_package.name}/data/tree95df.csv`);
    const view: DG.TableView = grok.shell.addTableView(df);
    await grok.functions.call(`${_package.name}:injectTreeToGrid`,
      {grid: view.grid, newickText: nwk, leafColName: 'id'});

    await awaitGrid(view.grid);
  });

  test('treeToGridApp', async () => {
    const app: TreeToGridApp = await grok.functions.call(`${_package.name}:treeToGridApp`, {});

    await awaitGrid(app.tableView!.grid);
  });

  test('treeCutAsTreeApp', async () => {
    await grok.functions.call(`${_package.name}:treeCutAsTreeApp`, {});
    await waitUntil(() => true, 'Error');
  });

  test('treeInGridCellApp', async () => {
    await grok.functions.call(`${_package.name}:treeInGridCellApp`, {});
    await waitUntil(() => true, 'Error');
  });

  test('getPhylocanvasGlService', async () => {
    const pglSvc: PhylocanvasGlService = await grok.functions.call(`${_package.name}:getPhylocanvasGlService`, {});
    await waitUntil(() => true, 'Error');
  });
});
