import * as DG from 'datagrok-api/dg';
import {Scope, signal} from '@datagrok-libraries/u2';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {BandChart, axisTitle, collectBands, lineSeries, niceCeil, niceFloor, niceStep, refLines, tooltip}
  from '../band-chart';
import {animalsByDay, animalsByDayFrame, lostLine} from '../analyze-view';

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

  test('A forecast axis tops at a round number', async () => {
    expect([1950, 2000, 2001, 1, 0.3, 74000].map(niceCeil).join(','), '2000,2000,2500,1,0.5,100000');
    // and its ticks divide that top (no 1800 then 2000)
    expect([2000, 2500, 1000, 5000, 0.5].map(niceStep).join(','), '500,500,200,1000,0.1');
    // a log axis' bottom: never the data's own 17.1209878087
    expect([17.1209878087, 10, 260, 0.37].map(niceFloor).join(','), '10,10,200,0.2');
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
    const owner = Scope.root('band chart test');
    const chart = owner.build(() => new BandChart(df, {x: 'day', lower: 'q05', upper: 'q95', line: 'median',
      split: 'arm', log}));
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
      owner.dispose();
    }
  });

  // Analyze's growth chart: the animals behind each point, and the stretch where fewer than 3 remain dashed
  const thinning = () => DG.DataFrame.fromCsv(
    `group,day,median,lower,upper,n
    Vehicle,0,100,90,110,10
    Vehicle,7,400,300,500,9
    Vehicle,14,900,700,1100,7
    Drug,0,100,90,110,10
    Drug,7,200,150,250,4
    Drug,14,260,200,300,2
    Drug,21,300,250,350,1`);

  test('animals by day: one line on who lost animals, the counts as a table', async () => {
    const a = animalsByDay(thinning(), ['Vehicle', 'Drug']);
    expect(lostLine(a), 'Animals by day: Vehicle lost 3, Drug lost 9');
    expect(lostLine(animalsByDay(DG.DataFrame.fromObjects([{group: 'A', day: 0, n: 5}, {group: 'A', day: 7, n: 5}])!,
      ['A'])), null);
    const t = animalsByDayFrame(a);
    expect(t.columns.names().join(), 'Group,Day 0,Day 7,Day 14,Day 21');
    expect(t.get('Day 14', 1), 2);
  });

  test('each point says how many animals it is made of; fewer than 3: its stretch dashed', async () => {
    const df = thinning();
    const [vehicle, drug] = collectBands(df, df.col('day')!, df.col('lower')!, df.col('upper')!, df.col('median'),
      df.col('group'), df.col('n'));
    expect(drug.n!.join(), '10,4,2,1');
    const fmt = {xUnits: 'd', yUnits: 'mm³', x: '', y: '', xTitle: 'Day'};
    const text = tooltip([{seriesName: 'Drug', dataIndex: 2, data: [14, 260]}], [vehicle, drug], fmt);
    expect(text.includes('260 mm³') && text.includes('2 of 10 animals'), true, text);
    // the vehicle never thins: one solid line
    const solid = lineSeries(vehicle, '#000', true, false) as any[];
    expect(solid.length, 1);
    expect(solid[0].lineStyle.type, 'solid');
    // the drug: solid to day 7, dashed from day 7 (4 animals) to day 21
    const s = lineSeries(drug, '#000', true, false) as any[];
    expect(s.length, 2);
    const shown = (x: any) => x.data.filter((p: any[]) => p[1] !== null).map((p: any[]) => p[0]).join();
    expect(shown(s[0]), '0,7');
    expect(s[1].lineStyle.type, 'dashed');
    expect(shown(s[1]), '7,14,21');
    // every series keeps the band's indices, for the tooltip
    expect(s.every((x) => x.data.length === 4), true);
  });
});
