import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, delay} from '@datagrok-libraries/test/src/test';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {getOrCreateParsedChartData} from '../fit/fit-chart-data';
import {multiSeriesCurveJson} from './curve-data';
import {MultiCurveViewer} from '../fit/multi-curve-viewer';

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
});
