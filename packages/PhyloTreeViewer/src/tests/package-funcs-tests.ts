import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package-test';
import {
  after,
  before,
  category,
  test,
  expect,
  expectArray,
  expectObject,
  awaitCheck
} from '@datagrok-libraries/utils/src/test';
import {PhylocanvasGlService} from '../utils/phylocanvas-gl-service';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TreeToGridApp} from '../apps/tree-to-grid-app';

/** Tests for package functions, test apps, file previews, file handlers, ... */
category('packageFuncs', () => {
  test('phylocanvasGlViewerApp', async () => {
    await grok.functions.call(`${_package.name}:phylocanvasGlViewerApp`, {});
    await awaitCheck(() => true, 'Error', 200);
  });

  test('injectTreeToGrid', async () => {
    const treeHelper: ITreeHelper = await getTreeHelper();
    const nwk = await grok.dapi.files.readAsText(`System:AppData/${_package.name}/data/tree95.nwk`);
    const df: DG.DataFrame = await grok.dapi.files.readCsv(`System.AppData/${_package.name}/data/tree95df.csv`);
    const view: DG.TableView = grok.shell.addTableView(df);
    await grok.functions.call(`${_package.name}:injectTreeToGrid`,
      {grid: view.grid, newickText: nwk, leafColName: 'id'});
    await awaitCheck(() => true, 'Error', 200);
  }, {skipReason: 'GROK-12877'});

  test('treeToGridApp', async () => {
    const app: TreeToGridApp = await grok.functions.call(`${_package.name}:treeToGridApp`, {});
    await awaitCheck(() => true, 'Error', 200);
  }, {skipReason: 'GROK-12877'});

  test('treeCutAsTreeApp', async () => {
    await grok.functions.call(`${_package.name}:treeCutAsTreeApp`, {});
    await awaitCheck(() => true, 'Error', 200);
  });

  test('treeInGridCellApp', async () => {
    await grok.functions.call(`${_package.name}:treeInGridCellApp`, {});
    await awaitCheck(() => true, 'Error', 200);
  });

  test('getPhylocanvasGlService', async () => {
    const pglSvc: PhylocanvasGlService = await grok.functions.call(`${_package.name}:getPhylocanvasGlService`, {});
    await awaitCheck(() => true, 'Error', 200);
  });
});
