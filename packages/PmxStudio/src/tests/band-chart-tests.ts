import * as DG from 'datagrok-api/dg';
import {signal} from '@datagrok-libraries/u2';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {BandChart, axisTitle, collectBands, refLines} from '../band-chart';

// The Plan and Analyze views' band chart: its bands by category and X, the table's reference lines, the axis titles,
// and the chart drawn over a table signal and its Log switch.
category('Band chart', () => {
  const table = () => DG.DataFrame.fromCsv(
    `day,arm,q05,median,q95
    2,control,110,150,210
    0,control,90,100,112
    1,control,98,120,150
    0,treated,90,100,112
    1,treated,95,110,130
    2,treated,99,118,150`);

  test('Bands by split, sorted by X', async () => {
    const df = table();
    const bands = collectBands(df, df.col('day')!, df.col('q05')!, df.col('q95')!, df.col('median'), df.col('arm'));
    expect(bands.map((b) => b.name).join(','), 'control,treated');
    expect(bands[0].x.join(','), '0,1,2');
    expect(bands[0].upper.join(','), '112,150,210');
    expect(bands[1].line!.join(','), '100,110,118');
  });

  test('Filter', async () => {
    const df = table();
    df.filter.set(5, false);
    const bands = collectBands(df, df.col('day')!, df.col('q05')!, df.col('q95')!, null, df.col('arm'));
    expect(bands[1].x.join(','), '0,1');
    expect(bands[1].line, null);
  });

  test('Reference lines', async () => {
    const lines = refLines([
      {formula: '${q95} = 1500', title: 'Humane limit'},
      {formula: '${day} = 19', title: 'Evaluation', style: 'solid'},
      {formula: '${q95} = ${median} * 2'},
      {formula: '${q05} = 100', visible: false},
    ], 'day');
    expect(lines.length, 2);
    expect(lines[0].vertical, false);
    expect(lines[0].value, 1500);
    expect(lines[0].dashed, true);
    expect(lines[1].vertical, true);
    expect(lines[1].dashed, false);
  });

  test('Axis titles carry the units once', async () => {
    expect(axisTitle('Tumor volume', 'mm³'), 'Tumor volume (mm³)');
    expect(axisTitle('Volume (mm³)', 'mm³'), 'Volume (mm³)');
    expect(axisTitle('day', ''), 'day');
    expect(axisTitle('Tumor volume', 'mm³', true), 'Tumor volume (mm³, log scale)');
  });

  test('Draws the table it is given, and none', async () => {
    const df = signal<DG.DataFrame | undefined>(undefined);
    const log = signal(true);
    const chart = new BandChart(df, {x: 'day', lower: 'q05', upper: 'q95', line: 'median', split: 'arm', log});
    document.body.append(chart.root);
    try {
      expect(chart.bands.length, 0);
      df.value = table();
      expect(chart.bands.map((b) => b.name).join(','), 'control,treated');
      log.value = false;
      expect(chart.bands.length, 2);
      df.value = undefined;
      expect(chart.bands.length, 0);
    } finally {
      chart.root.remove();
      chart.dispose();
    }
  });
});
