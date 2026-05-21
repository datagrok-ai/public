import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, look} from '../helpers';

// Regression coverage for GROK-19466 (1.27.0): Trellis Bar chart with
// globalScale=true used to break when valueAggrType was swapped under the
// inner BarChartLook. The fix preserves viewer state across all aggr swaps.
//
// Detached viewers throughout (no TableView, no canvas geometry); state is
// pinned via getOptions(true).look. innerViewerLook is typed `any` and
// backed by a Dart LookAndFeel rebuilt on every reassignment, so we read
// back via the serialized JSON envelope rather than via props.
const AGGRS = ['avg', 'sum', 'min', 'max', 'count', 'med', 'q1', 'q3', 'stdev'];

function innerLook(aggr: string): {[k: string]: any} {
  return {'#type': 'BarChartLook', valueColumnName: 'age', valueAggrType: aggr, splitColumnName: 'race'};
}

function makeTrellis(): DG.Viewer {
  return DG.Viewer.trellisPlot(demog(40), {
    viewerType: DG.VIEWER.BAR_CHART, xColumnNames: ['race'], yColumnNames: ['sex'],
    globalScale: true, innerViewerLook: innerLook('avg'),
  });
}

category('AI: GROK-19466: Trellis Bar chart global scale + aggr', () => {
  test('construct trellis Bar chart with globalScale=true', async () => {
    var v: DG.Viewer | null = null;
    expectNoThrow(() => { v = makeTrellis(); });
    expect(v!.type, DG.VIEWER.TRELLIS_PLOT);
    expectLook(v!, {viewerType: DG.VIEWER.BAR_CHART, globalScale: true});
    expect(look(v!)['innerViewerLook']['valueAggrType'], 'avg');
    expect(look(v!)['innerViewerLook']['valueColumnName'], 'age');
  });

  test('switch valueAggrType to avg under globalScale', async () => {
    const v = makeTrellis();
    expectNoThrow(() => v.setOptions({innerViewerLook: innerLook('avg')}));
    expect(look(v)['innerViewerLook']['valueAggrType'], 'avg');
  });

  test('cycle through every aggr type under globalScale=true', async () => {
    const v = makeTrellis();
    for (var aggr of AGGRS) {
      expectNoThrow(() => v.setOptions({innerViewerLook: innerLook(aggr)}));
      expect(look(v)['innerViewerLook']['valueAggrType'], aggr);
    }
  });

  test('toggle globalScale between true/false after each aggr change', async () => {
    const v = makeTrellis();
    var scale = true;
    for (var aggr of AGGRS) {
      expectNoThrow(() => {
        v.setOptions({innerViewerLook: innerLook(aggr)});
        scale = !scale; v.setOptions({globalScale: scale});
        scale = !scale; v.setOptions({globalScale: scale});
      });
      expectLook(v, {globalScale: scale});
      expect(look(v)['innerViewerLook']['valueAggrType'], aggr);
    }
  });

  test('globalScale=false then switch aggr then globalScale=true', async () => {
    const v = makeTrellis();
    expectNoThrow(() => {
      v.setOptions({globalScale: false});
      v.setOptions({innerViewerLook: innerLook('sum')});
      v.setOptions({globalScale: true});
      v.setOptions({innerViewerLook: innerLook('max')});
    });
    expectLook(v, {globalScale: true});
    expect(look(v)['innerViewerLook']['valueAggrType'], 'max');
  });
});
