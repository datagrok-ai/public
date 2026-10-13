import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray, expectObject} from '@datagrok-libraries/test/src/test';
import {chartSpec, ChartSpec} from '../components/RunComparison/chart-spec';
import {Comparison} from '../components/RunComparison/comparison-builders';

// the spec reads only the row count of the frame
const frame = (rowCount = 3) => DG.DataFrame.create(rowCount);

const scalar = (rowCount = 3): Comparison =>
  ({kind: 'scalar', chartDf: frame(rowCount), valueColumnName: 'X (g/L)'});

const multiScalar = (valueColumnNames: string[]): Comparison =>
  ({kind: 'multi-scalar', chartDf: frame(), valueColumnNames});

const column = (fields: {
  valueColumnNames?: string[], splitColumnName?: string, isKeyIndex?: boolean, isScatter?: boolean,
  melted?: {seriesColumnName: string, valueColumnName: string},
} = {}): Comparison => {
  const valueColumnNames = fields.valueColumnNames ?? ['Y'];
  return {
    kind: 'column', chartDf: frame(), valueColumnName: valueColumnNames[0], indexColumnName: 'Time (h)',
    isKeyIndex: false, isScatter: false, ...fields, valueColumnNames,
  };
};

const viewerSpec = (spec: ChartSpec) => {
  if (spec.kind !== 'viewer')
    throw new Error(`expected a viewer, got the hint "${spec.text}"`);
  return spec;
};

const hintText = (spec: ChartSpec) => spec.kind === 'hint' ? spec.text : null;

category('RunComparison: chart spec hints', () => {
  test('no selection asks to select a value', async () => {
    expect(hintText(chartSpec(scalar(), false, 'radar', true)), 'Select what to compare above');
  });

  test('a missing or empty comparison has nothing to show', async () => {
    const nothing = 'Nothing to show: no data points across the selected runs';
    expect(hintText(chartSpec(null, true, 'radar', true)), nothing);
    expect(hintText(chartSpec(scalar(0), true, 'radar', true)), nothing);
  });

  test('multiple values over a key index ask for another index', async () => {
    const text = hintText(chartSpec(column({isKeyIndex: true, valueColumnNames: ['Y', 'Z']}), true, 'radar', true));
    expect(text?.startsWith('Multiple values cannot chart over a categorical index'), true);
  });
});

category('RunComparison: chart spec scalars', () => {
  test('a scalar charts as a bar per run', async () => {
    const spec = viewerSpec(chartSpec(scalar(), true, 'radar', true));
    expect(spec.type, DG.VIEWER.BAR_CHART);
    expectObject(spec.options, {valueColumnName: 'X (g/L)', valueAggrType: 'avg', splitColumnName: 'Run'});
    expect(spec.minHeight, 250);
    expect(spec.scalarChartSwitch, false);
  });

  test('two scalars chart as a PC plot without the switch', async () => {
    const spec = viewerSpec(chartSpec(multiScalar(['A', 'B']), true, 'radar', true));
    expect(spec.type, DG.VIEWER.PC_PLOT);
    expectArray(spec.options.columnNames, ['A', 'B']);
    expect(spec.options.colorColumnName, 'Run');
    expect(spec.minHeight, 300);
    expect(spec.scalarChartSwitch, false);
  });

  test('three scalars chart as a radar with the switch', async () => {
    const spec = viewerSpec(chartSpec(multiScalar(['A', 'B', 'C']), true, 'radar', true));
    expect(spec.type, 'Radar');
    expectArray(spec.options.valuesColumnNames, ['A', 'B', 'C']);
    expect(spec.options.colorColumnName, 'Run');
    expect(spec.options.showCurrentRow, false);
    expect(spec.minHeight, 300);
    expect(spec.scalarChartSwitch, true);
  });

  test('a picked PC plot keeps the switch', async () => {
    const spec = viewerSpec(chartSpec(multiScalar(['A', 'B', 'C']), true, 'pcplot', true));
    expect(spec.type, DG.VIEWER.PC_PLOT);
    expect(spec.scalarChartSwitch, true);
  });

  test('without the radar viewer three scalars fall back to a PC plot', async () => {
    const spec = viewerSpec(chartSpec(multiScalar(['A', 'B', 'C']), true, 'radar', false));
    expect(spec.type, DG.VIEWER.PC_PLOT);
    expect(spec.scalarChartSwitch, false);
  });
});

category('RunComparison: chart spec columns', () => {
  test('a key index charts as bars split by index and stacked by run', async () => {
    const spec = viewerSpec(chartSpec(column({isKeyIndex: true}), true, 'radar', true));
    expect(spec.type, DG.VIEWER.BAR_CHART);
    expectObject(spec.options,
      {valueColumnName: 'Y', valueAggrType: 'avg', splitColumnName: 'Time (h)', stackColumnName: 'Run'});
    expect(spec.minHeight, 250);
  });

  test('a single value charts as a line split by run', async () => {
    const spec = viewerSpec(chartSpec(column(), true, 'radar', true));
    expect(spec.type, DG.VIEWER.LINE_CHART);
    expect(spec.options.xColumnName, 'Time (h)');
    expectArray(spec.options.yColumnNames, ['Y']);
    expectArray(spec.options.splitColumnNames, ['Run']);
    expect(spec.options.multiAxis, false);
    expect(spec.minHeight, 250);
  });

  test('a split column adds an inner line split', async () => {
    const spec = viewerSpec(chartSpec(column({splitColumnName: 'Reactor'}), true, 'radar', true));
    expectArray(spec.options.splitColumnNames, ['Run', 'Reactor']);
  });

  test('multiple values share one line chart, taller per value', async () => {
    const spec = viewerSpec(chartSpec(column({valueColumnNames: ['Y', 'Z', 'W']}), true, 'radar', true));
    expect(spec.type, DG.VIEWER.LINE_CHART);
    expectArray(spec.options.yColumnNames, ['Y', 'Z', 'W']);
    expect(spec.minHeight, 750);
  });

  test('independent points color by run', async () => {
    const spec = viewerSpec(chartSpec(column({isScatter: true}), true, 'radar', true));
    expect(spec.type, DG.VIEWER.SCATTER_PLOT);
    expectObject(spec.options, {xColumnName: 'Time (h)', yColumnName: 'Y', colorColumnName: 'Run'});
    expect('markersColumnName' in spec.options, false);
    expect(spec.minHeight, 300);
  });

  test('a split takes the color and runs move to markers', async () => {
    const spec = viewerSpec(chartSpec(column({isScatter: true, splitColumnName: 'Reactor'}), true, 'radar', true));
    expectObject(spec.options, {colorColumnName: 'Reactor', markersColumnName: 'Run'});
  });

  test('melted values color by value name with runs as markers', async () => {
    const spec = viewerSpec(chartSpec(column({
      isScatter: true, splitColumnName: 'Reactor', valueColumnNames: ['Y', 'Z'],
      melted: {seriesColumnName: 'Data', valueColumnName: 'Value'},
    }), true, 'radar', true));
    expect(spec.type, DG.VIEWER.SCATTER_PLOT);
    expectObject(spec.options,
      {xColumnName: 'Time (h)', yColumnName: 'Value', colorColumnName: 'Data', markersColumnName: 'Run'});
    expect(spec.minHeight, 300);
  });
});
