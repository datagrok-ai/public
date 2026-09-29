import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {awaitCheck, category, delay, test, testViewer} from '@datagrok-libraries/test/src/test';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {NodeType} from '@datagrok-libraries/bio/src/trees';
import {parseNewick} from '@datagrok-libraries/bio/src/trees/phylocanvas';

import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';
import {TreeHelper} from '../utils/tree-helper';

import {_package} from '../package-test';

category('GridWithTree', () => {
  test('open', async () => {
    const _th: ITreeHelper = new TreeHelper();

    const csv: string = await _package.files.readAsText('data/tree95df.csv');
    const newickStr: string = await _package.files.readAsText('data/tree95.nwk');
    const leafColName = 'id';

    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const newickRoot: NodeType = parseNewick(newickStr);

    const tv: DG.TableView = grok.shell.addTableView(dataDf);
    await awaitCheck(() => {
      return $(tv.root).find('.d4-grid canvas').length > 0;
    }, 'The view grid canvas not found', 200);
    ;const neighborWidth = 250;
    injectTreeForGridUI2(tv.grid, newickRoot, leafColName, neighborWidth);
    await awaitCheck(() => {
      return $(tv.root).find('.ui-div canvas').length == 1;
    }, 'Injected tree not found', 200);
  });

  test('closedView', async () => {
    const csv: string = await _package.files.readAsText('data/tree95df.csv');
    const newickStr: string = await _package.files.readAsText('data/tree95.nwk');
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);

    const tv: DG.TableView = grok.shell.addTableView(dataDf);
    await awaitCheck(() => {
      return $(tv.root).find('.d4-grid canvas').length > 0;
    }, 'The view grid canvas not found', 200);
    const grid: DG.Grid = tv.grid;
    injectTreeForGridUI2(grid, parseNewick(newickStr), 'id', 250);

    tv.close();
    await awaitCheck(() => grid.dataFrame == null, 'The grid is still attached to the table', 1000);

    grok.shell.clearLastError();
    dataDf.mouseOverRowIdx = 1;
    dataDf.currentRowIdx = 2;
    await delay(500);
    const unhandled = await grok.shell.lastError;
    if (unhandled)
      throw new Error(unhandled);
  });
});
