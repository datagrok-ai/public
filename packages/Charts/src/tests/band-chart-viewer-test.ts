import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getOptions} from './utils';
import {axisTitle, collectBands, refLines} from '../viewers/band-chart/band-chart-viewer';


category('Band chart', () => {
  const TYPE = 'Band chart';
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = DG.DataFrame.fromCsv(
      `day, arm, q05, median, q95
      2, control, 110, 150, 210
      0, control, 90, 100, 112
      1, control, 98, 120, 150
      0, treated, 90, 100, 112
      1, treated, 95, 110, 130
      2, treated, 99, 118, 150`);
    tv = grok.shell.addTableView(df);
  });

  test('Creation', async () => {
    const viewer = DG.Viewer.fromType(TYPE, df);
    expect(viewer instanceof DG.Viewer, true);
    expect(viewer.type, TYPE);
    expect(viewer.table.id, df.id);
  });

  test('Standard Properties', async () => {
    const viewer = DG.Viewer.fromType(TYPE, df);
    const options = await getOptions(viewer);
    // no band declared: the first three numerical columns
    expect(options.xColumnName, 'day');
    expect(options.lowerColumnName, 'q05');
    expect(options.upperColumnName, 'median');
    expect(options.yAxisType, 'linear');
    expect(options.bandOpacity, 0.25);
  });

  test('Bands by split, sorted by X', async () => {
    const bands = collectBands(df, df.col('day')!, df.col('q05')!, df.col('q95')!, df.col('median'), df.col('arm'));
    expect(bands.map((b) => b.name).join(','), 'control,treated');
    expect(bands[0].x.join(','), '0,1,2');
    expect(bands[0].upper.join(','), '112,150,210');
    expect(bands[1].line!.join(','), '100,110,118');
  });

  test('Filter', async () => {
    df.filter.set(5, false);
    const bands = collectBands(df, df.col('day')!, df.col('q05')!, df.col('q95')!, null, df.col('arm'));
    df.filter.setAll(true);
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
    expect(axisTitle('Tumour volume', 'mm³'), 'Tumour volume (mm³)');
    expect(axisTitle('Volume (mm³)', 'mm³'), 'Volume (mm³)');
    expect(axisTitle('day', ''), 'day');
    expect(axisTitle('Tumour volume', 'mm³', true), 'Tumour volume (mm³, log scale)');
    expect(axisTitle('Volume', '', true), 'Volume (log scale)');
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
