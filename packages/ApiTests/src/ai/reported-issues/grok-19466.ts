import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

const AGGRS = ['avg', 'sum', 'min', 'max', 'count', 'med', 'q1', 'q3', 'stdev'];

function makeTrellis(): DG.Viewer {
  const df = grok.data.demo.demog(40);
  return DG.Viewer.trellisPlot(df, {
    viewerType: DG.VIEWER.BAR_CHART,
    xColumnNames: ['race'],
    yColumnNames: ['sex'],
    globalScale: true,
    innerViewerLook: {
      '#type': 'BarChartLook',
      valueColumnName: 'age',
      valueAggrType: 'avg',
      splitColumnName: 'race',
    },
  });
}

category('AI: GROK-19466: Trellis Bar chart global scale + aggr', () => {
  test('construct trellis Bar chart with globalScale=true', async () => {
    var threw = false;
    var v: DG.Viewer | null = null;
    try {
      v = makeTrellis();
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v != null, true);
    expect(v!.type, DG.VIEWER.TRELLIS_PLOT);
    const look = v!.getOptions(true).look;
    expect(look['viewerType'], DG.VIEWER.BAR_CHART);
    expect(look['globalScale'], true);
    expect(look['innerViewerLook'] != null, true);
    expect(look['innerViewerLook']['valueAggrType'], 'avg');
    expect(look['innerViewerLook']['valueColumnName'], 'age');
  });

  test('switch valueAggrType to avg under globalScale', async () => {
    const v = makeTrellis();
    var threw = false;
    try {
      v.setOptions({innerViewerLook: {
        '#type': 'BarChartLook',
        valueColumnName: 'age',
        valueAggrType: 'avg',
        splitColumnName: 'race',
      }});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    const inner = v.getOptions(true).look['innerViewerLook'];
    expect(inner != null, true);
    expect(inner['valueAggrType'], 'avg');
  });

  test('cycle through every aggr type under globalScale=true', async () => {
    const v = makeTrellis();
    for (var aggr of AGGRS) {
      var threw = false;
      try {
        v.setOptions({innerViewerLook: {
          '#type': 'BarChartLook',
          valueColumnName: 'age',
          valueAggrType: aggr,
          splitColumnName: 'race',
        }});
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      const inner = v.getOptions(true).look['innerViewerLook'];
      expect(inner != null, true);
      expect(inner['valueAggrType'], aggr);
    }
  });

  test('toggle globalScale between true/false after each aggr change', async () => {
    const v = makeTrellis();
    var scale = true;
    for (var aggr of AGGRS) {
      var threw = false;
      try {
        v.setOptions({innerViewerLook: {
          '#type': 'BarChartLook',
          valueColumnName: 'age',
          valueAggrType: aggr,
          splitColumnName: 'race',
        }});
        scale = !scale;
        v.setOptions({globalScale: scale});
        scale = !scale;
        v.setOptions({globalScale: scale});
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      const look = v.getOptions(true).look;
      expect(look['globalScale'], scale);
      expect(look['innerViewerLook']['valueAggrType'], aggr);
    }
  });

  test('globalScale=false then switch aggr then globalScale=true', async () => {
    const v = makeTrellis();
    var threw = false;
    try {
      v.setOptions({globalScale: false});
      v.setOptions({innerViewerLook: {
        '#type': 'BarChartLook',
        valueColumnName: 'age',
        valueAggrType: 'sum',
        splitColumnName: 'race',
      }});
      v.setOptions({globalScale: true});
      v.setOptions({innerViewerLook: {
        '#type': 'BarChartLook',
        valueColumnName: 'age',
        valueAggrType: 'max',
        splitColumnName: 'race',
      }});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    const look = v.getOptions(true).look;
    expect(look['globalScale'], true);
    expect(look['innerViewerLook']['valueAggrType'], 'max');
  });
});
