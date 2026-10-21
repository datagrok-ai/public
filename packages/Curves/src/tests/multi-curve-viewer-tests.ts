import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {category, test, expect, delay, testEvent} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {getOrCreateParsedChartData} from '../fit/fit-chart-data';
import {MultiCurveViewer} from '../fit/multi-curve-viewer';
import {multiSeriesCurveJson} from './curve-data';

function curveTable(name: string): DG.DataFrame {
  const col = DG.Column.fromStrings('curve', [multiSeriesCurveJson([-7, -5]), multiSeriesCurveJson([-6, -4])]);
  col.semType = FitConstants.FIT_SEM_TYPE;
  const df = DG.DataFrame.fromColumns([col]);
  df.name = name;
  return df;
}

category('multi curve viewer', () => {
  test('merging cell series leaves the cell alone, and unmerging comes back', async () => {
    const df = curveTable('multiCurveMerge');
    const tv = grok.shell.addTableView(df);
    try {
      df.currentRowIdx = 0;
      tv.addViewer('MultiCurveViewer', {curvesColumnNames: ['curve']});
      await delay(1000);
      // addViewer hands back the host wrapper, and the merge happens on the JS instance behind it
      const viewer = Array.from(tv.viewers).find((v) => v.type === 'MultiCurveViewer') as MultiCurveViewer;
      expect(viewer !== undefined, true, 'the viewer was not added');
      expect(viewer.data?.series?.length, 2, 'the current row carries two curves');

      viewer.props.set('mergeCellSeries', true as unknown as object);
      await delay(500);
      expect(viewer.data?.series?.length, 1, 'merging should leave one curve in the viewer');
      // the parsed data is cached and shared with the grid - merging into it used to rewrite the cell
      expect(getOrCreateParsedChartData(df.cell(0, 'curve')).series!.length, 2,
        'the cell itself must still hold both curves');

      viewer.props.set('mergeCellSeries', false as unknown as object);
      await delay(500);
      expect(viewer.data?.series?.length, 2, 'unmerging should bring both curves back');
    } finally {
      tv.close();
    }
  }, {timeout: 60000});

  test('the viewer names where each of its curves came from', async () => {
    const df = curveTable('multiCurveFitIdentities');
    const tv = grok.shell.addTableView(df);
    try {
      df.currentRowIdx = 0;
      tv.addViewer('MultiCurveViewer', {curvesColumnNames: ['curve']});
      await delay(1500);
      const viewer = Array.from(tv.viewers).find((v) => v.type === 'MultiCurveViewer') as any;
      expect(viewer.fitIdentities.length, viewer.data.series.length,
        'every series should carry an identity, in the order they are drawn');
      expect(viewer.fitIdentities.every((id: string | undefined) => id !== undefined), true,
        'and none of them should be missing');
    } finally {
      tv.close();
    }
  }, {timeout: 60000});

  test('the viewer answers the same hovers the grid does', async () => {
    const logIC50s = Array.from({length: 40}, (_, i) => -8 + i * 0.1);
    const col = DG.Column.fromStrings('curve', [multiSeriesCurveJson(logIC50s)]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'multiCurveHover';
    const tv = grok.shell.addTableView(df);
    try {
      df.currentRowIdx = 0;
      tv.addViewer('MultiCurveViewer', {curvesColumnNames: ['curve']});
      await delay(2000);
      const viewer = Array.from(tv.viewers).find((v) => v.type === 'MultiCurveViewer') as any;
      expect(viewer.data.series.length, 40, 'the current row carries all the curves');

      // the grid routes mouse moves to the cell renderer; the viewer has to do it itself
      const canvas: HTMLCanvasElement = viewer.canvas;
      const rect = canvas.getBoundingClientRect();
      // any row will do - whether the tail collapses into `+N more` depends on the viewer's height
      let tooltip: string[] = [];
      for (let x = 10; x < canvas.width - 5 && tooltip.length === 0; x += 10) {
        for (let y = 6; y < canvas.height - 5; y += 6) {
          canvas.dispatchEvent(new MouseEvent('mousemove',
            {bubbles: true, clientX: rect.left + x, clientY: rect.top + y}));
          const shown = document.querySelector('.d4-tooltip') as HTMLElement | null;
          const lines: string[] = shown?.innerText ?
            shown.innerText.split(String.fromCharCode(10)).filter((t: string) => t.trim().length > 0) : [];
          if (lines.some((line) => line.startsWith('series '))) {
            tooltip = lines;
            break;
          }
        }
      }
      expect(tooltip.length > 0, true, 'hovering the legend in the viewer should name a curve');
    } finally {
      tv.close();
    }
  }, {timeout: 90000});

  // GROK-20661: a range one cell declares for its own chart was enforced on the whole overlay, cutting the rest
  test('the axes take in every point, whatever range one cell declares for its own chart', async () => {
    const curve = (name: string, value: object) => {
      const col = DG.Column.fromStrings(name, [JSON.stringify(value)]);
      col.semType = FitConstants.FIT_SEM_TYPE;
      return col;
    };
    const df = DG.DataFrame.fromColumns([
      curve('narrow', {chartOptions: {minY: 0, maxY: 10}, series: [{points: [{x: 0.1, y: 1}, {x: 1, y: 5}, {x: 10, y: 9}]}]}),
      curve('tall', {series: [{points: [{x: 0.1, y: 10}, {x: 1, y: 50}, {x: 10, y: 100}]}]}),
    ]);
    df.name = 'multiCurveAxes';
    const tv = grok.shell.addTableView(df);
    try {
      df.currentRowIdx = 0;
      tv.addViewer('MultiCurveViewer', {curvesColumnNames: ['narrow', 'tall']});
      await delay(2000);
      const viewer = Array.from(tv.viewers).find((v) => v.type === 'MultiCurveViewer') as MultiCurveViewer;
      const values = viewer.getWidgetStatus().values!;
      expect(values['curves'], 2);
      expect(values['points outside the axes'], 0, `the axes span ${values['y range']}`);
    } finally {
      tv.close();
    }
  });

  test('the viewer reports the legend it drew, and a render still on its way', async () => {
    const col = DG.Column.fromStrings('curve', [multiSeriesCurveJson(Array.from({length: 40}, (_, i) => -8 + i * 0.1))]);
    col.semType = FitConstants.FIT_SEM_TYPE;
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'multiCurveStatus';
    const tv = grok.shell.addTableView(df);
    try {
      df.currentRowIdx = 0;
      tv.addViewer('MultiCurveViewer', {curvesColumnNames: ['curve']});
      await delay(2000);
      const viewer = Array.from(tv.viewers).find((v) => v.type === 'MultiCurveViewer') as MultiCurveViewer;
      let status = viewer.getWidgetStatus();
      expect(status.parts['canvas'] === viewer.canvas, true, 'the canvas is the part the hit areas are on');
      expect(status.values!['curves'], 40);
      expect(status.values!['legend shown'], true);
      expect(viewer.props.get('showLegend') as unknown, true, 'Show Legend should read on while the legend is drawn');
      expect((status.values!['legend rows'] as number) + (status.values!['legend hidden'] as number), 40,
        'the rows drawn and the ones counted should cover every curve');
      const more = status.hitAreas['legend more'];
      expect(more !== undefined, true, 'forty curves do not fit, so the rest is counted');

      // the area is where the hover answers with the list
      const rect = viewer.canvas.getBoundingClientRect();
      viewer.canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true,
        clientX: rect.left + more.x + more.width / 2, clientY: rect.top + more.y + more.height / 2}));
      const listed = document.querySelector('.d4-tooltip')?.querySelectorAll('canvas').length ?? 0;
      ui.tooltip.hide();
      expect(listed > 1, true, `hovering "legend more" should list the curves, it listed ${listed}`);

      viewer.props.set('showLegend', false as unknown as object);
      await delay(500);
      status = viewer.getWidgetStatus();
      expect(status.values!['legend shown'], false, 'Show Legend off should leave the curves unnamed');
      expect(Object.keys(status.hitAreas).some((name) => name.includes('legend')), false,
        'and no legend row to point at');

      let pendingOnChange = false;
      const sub = df.onSelectionChanged.subscribe(() => pendingOnChange = viewer.isRenderPending);
      try {
        await testEvent(viewer.onRendered, () => expect(viewer.isRenderPending, false, 'a render landed'),
          () => df.selection.set(0, true), 5000);
      } finally {
        sub.unsubscribe();
      }
      expect(pendingOnChange, true, 'a selection change should say a render is on its way');
    } finally {
      tv.close();
    }
  }, {timeout: 90000});
});
