import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {fromEvent} from 'rxjs';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';

category('rendering-tooltip', () => {
  test('scatterPlot', async () => {
    await _testScatterPlotTooltip();
  });

  const molCoordsCsv = `mol,x,y
[H:1]C(=O)O,0,0
[H:2]N[H:1],1,0
[H:1]NC(C#C)C([OH:2])=O,0,2
[H:1]NC(CN=[N+]=[N-])C([OH:2])=O,1,1`;

  async function _testScatterPlotTooltip(): Promise<void> {
    const df = DG.DataFrame.fromCsv(molCoordsCsv);
    df.currentRowIdx = 0;
    const view = grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);
    const sp: DG.ScatterPlotViewer = df.plot.scatter({x: 'x', y: 'y'});
    view.dockManager.dock(sp, DG.DOCK_TYPE.RIGHT, null);
    await delay(200);
    await Promise.all([
      testEvent(sp.onAfterDrawScene, () => {}, () => { sp.invalidateCanvas(); }, 1000),
      testEvent(view.grid.onAfterDrawContent, () => {}, () => { view.grid.invalidate(); }, 1000)
    ]);

    const spBcr = sp.root.getBoundingClientRect();
    const wp = sp.worldToScreen(1, 0);
    const ev = new MouseEvent('mousemove', {
      cancelable: true, bubbles: true, view: window, button: 0,
      clientX: spBcr.left + wp.x, clientY: spBcr.top + wp.y
    });

    const spCanvas = $(sp.root).find('canvas').get()[0] as HTMLCanvasElement;
    await testEvent(fromEvent(spCanvas, 'mousemove'), () => {
      _package.logger.debug(`Test: event, currentRowIdx: ${df.currentRowIdx}`);
      _package.logger.debug(`Test: event, tooltip.root:\n${ui.tooltip.root.outerHTML}`);
      expect($(ui.tooltip.root).find('div table.d4-row-tooltip-table tr td canvas').length, 1);
      expect(sp.hitTest(wp.x, wp.y), 1);
    }, () => {
      spCanvas.dispatchEvent(ev);
    }, 500);
    // TODO: Any error occurred become 'Cannot read properties of null (reading 'get$columns')' because of scatter plot
    testEvent(view.grid.onAfterDrawContent, () => {}, () => { view.grid.invalidate(); }, 1000);
  }
});
