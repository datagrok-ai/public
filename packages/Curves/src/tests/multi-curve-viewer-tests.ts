import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, delay} from '@datagrok-libraries/test/src/test';
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
});
