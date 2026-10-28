import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {addSavedAnalyses, hideJmpHiddenColumns, jmpTableToDataFrame, savedAnalyses, skippedColumns}
  from '../formats/jmp/jmp-to-dataframe';
import {readJmp} from '../formats/jmp/jmp-parser';
import {JslNode, parseJsl} from '../formats/jmp/jsl';
import {fixture, fixtureBytes, FIXTURES} from './jmp-reader-tests';
import {makeLargeJmp, rowsForSize} from './large-jmp';

type JslCall = Extract<JslNode, {type: 'call'}>;

category('JMP: dataframe', () => {
  test('types and nulls', async () => {
    const df = jmpTableToDataFrame(await fixture('example1'), 'example1');
    expect(df.rowCount, 4);
    expect(df.col('ints')!.type, DG.COLUMN_TYPE.INT);
    expect(df.col('floats')!.type, DG.COLUMN_TYPE.FLOAT);
    expect(df.col('charconstwidth')!.type, DG.COLUMN_TYPE.STRING);
    // JMP text that reads as numbers stays text
    expect(df.col('formula')!.type, DG.COLUMN_TYPE.STRING);
    expectArray(df.col('formula')!.toList(), ['2', '4', '6', '8']);
    expect(df.col('time')!.type, DG.COLUMN_TYPE.DATE_TIME);
    expect(df.col('date')!.type, DG.COLUMN_TYPE.DATE_TIME);
    expect(df.col('date')!.meta.format, 'yyyy-MM-dd');
    expect(df.col('duration')!.type, DG.COLUMN_TYPE.FLOAT);
    expect(df.col('duration')!.meta.units, 's');
    expect(df.col('pressures')!.isNone(1), true);
    expect(df.col('time')!.isNone(3), true);
    expect(df.col('pressures')!.get(3), 4.63309e110);
    expect(df.col('char utf8')!.get(1), '🚴💨');
    expect(df.col('time')!.get(0).valueOf(), Date.UTC(1976, 3, 1, 21, 12));
  });

  test('integer nulls', async () => {
    const df = jmpTableToDataFrame(await fixture('byteintegers_notcompressed'));
    const four = df.col('fourbyte')!;
    expect(four.type, DG.COLUMN_TYPE.INT);
    expect(four.isNone(0), true);
    expect(four.get(1), 2147483647);
    expect(four.get(2), -2147483646);
  });

  test('table tags', async () => {
    const df = jmpTableToDataFrame(await fixture('example1'));
    expect(df.getTag('jmp.version'), '17.2.0');
    expect(df.getTag('jmp.saveTime'), '2024-03-02T08:06:26.000Z');
    expect(df.col('floats')!.getTag('jmp.format'), '01 00 0c 63 63 08');
  });

  test('geographic semantic types', async () => {
    const df = jmpTableToDataFrame(await fixture('geographic'));
    expect(df.col('Latitude_DDD')!.semType, DG.SEMTYPE.LATITUDE);
    expect(df.col('Longitude_DDD')!.semType, DG.SEMTYPE.LONGITUDE);
  });

  test('row states: a marker and a color column each', async () => {
    const t = await fixture('rowstate');
    const df = jmpTableToDataFrame(t);
    expect(df.rowCount, 3);
    expect(skippedColumns(t).length, 0);
    expectArray(df.columns.names(), ['rowstate1 marker', 'rowstate1 color', 'rowstate2 marker', 'rowstate2 color',
      'rowstate3 marker', 'rowstate3 color']);
    expectArray(df.col('rowstate2 marker')!.toList(), ['X', 'Y', 'Z']);
    expectArray(df.col('rowstate2 color')!.toList(), ['#000000', '#000000', '#C0C0C0']);
    expect(df.col('rowstate2 color')!.meta.colors.getType(), DG.COLOR_CODING_TYPE.CATEGORICAL);
    expect(JSON.parse(df.col('rowstate2 color')!.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL))['#C0C0C0'], '#C0C0C0');
    expect(JSON.parse(df.col('rowstate2 marker')!.getTag(DG.TAGS.MARKER_CODING))['X'], 'cross x border');
    expect(JSON.parse(df.col('rowstate3 marker')!.getTag(DG.TAGS.MARKER_CODING))['▲'], 'triangle top');
    expect(df.col('rowstate1 marker')!.getTag('jmp.rowState'), 'rowstate1');
  });

  test('saved analyses open as viewers', async () => {
    const t = await readJmp(await _package.files.readAsBytes('tests/saved-bar-chart.jmp'));
    const df = jmpTableToDataFrame(t, 'saved-bar-chart');
    expect(JSON.parse(df.getTag('jmp.scripts')!)[0].name, 'Bar chart of A');
    const tv = grok.shell.addTableView(df);
    try {
      expect(addSavedAnalyses(tv, savedAnalyses(t)), 1);
      const bar = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.BAR_CHART);
      expect(bar !== undefined, true);
      expect(bar!.props.splitColumnName, 'A');
    } finally {
      tv.close();
    }
  });

  test('a fit per Overlay group, from a JMP 15 file', async () => {
    const t = await readJmp(await _package.files.readAsBytes('tests/jmp15-graph-builder.jmp'));
    const tv = grok.shell.addTableView(jmpTableToDataFrame(t, 'jmp15-graph-builder'));
    try {
      expect(addSavedAnalyses(tv, savedAnalyses(t)), 1);
      const sp = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT)!;
      expect(sp.props.colorColumnName, 'Condition');
      expect(sp.props.showRegressionLine, true);
      expect(sp.props.regressionPerCategory, true);
    } finally {
      tv.close();
    }
  });

  test('a trellis of a converted Graph Builder', async () => {
    const df = grok.data.demo.demog(100);
    const tv = grok.shell.addTableView(df);
    try {
      const added = addSavedAnalyses(tv, [{name: 'Panels', text: '', notes: [], viewers: [{type: DG.VIEWER.TRELLIS_PLOT,
        options: {viewerType: DG.VIEWER.SCATTER_PLOT, xColumnNames: ['sex'], yColumnNames: [],
          innerViewerLook: {xColumnName: 'height', yColumnName: 'weight'}, title: 'Panels'}}]}]);
      expect(added, 1);
      const trellis = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.TRELLIS_PLOT)!;
      expect(trellis.props.viewerType, DG.VIEWER.SCATTER_PLOT);
      expectArray(trellis.props.xColumnNames, ['sex']);
      const look = JSON.stringify(trellis.getOptions(true));
      expect(look.includes('"xColumnName":"height"'), true, look.slice(0, 400));
    } finally {
      tv.close();
    }
  });

  test('currency codes and hidden columns', async () => {
    const cur = jmpTableToDataFrame(await fixture('currencies'));
    expect(cur.col('EUR')!.meta.units, 'EUR');
    expect(cur.col('USD')!.getTag('jmp.currency'), 'USD');
    const ints = jmpTableToDataFrame(await fixture('byteintegers'));
    expect(ints.col('Column 1')!.getTag('jmp.hidden'), 'true');
    expect(ints.col('2-byte integer')!.getTag('jmp.hidden'), null);
    const tv = grok.shell.addTableView(ints);
    try {
      hideJmpHiddenColumns(tv);
      expect(tv.grid.columns.byName('Column 1')!.visible, false);
      expect(tv.grid.columns.byName('2-byte integer')!.visible, true);
    } finally {
      tv.close();
    }
  });

  test('column properties: units, formula, notes, value order', async () => {
    const t = await fixture('example1');
    const chars = t.columns.find((c) => c.info.name === 'charconstwidth')!;
    chars.properties = {properties: ['Notes("Letters")', 'Value Order({Custom Order({"d", "c", "b", "a"})})']
      .map((text) => ({name: text.slice(0, text.indexOf('(')), text, node: parseJsl(text) as JslCall}))};
    const df = jmpTableToDataFrame(t);
    expect(df.col('pressures')!.meta.units, 'kPa');
    expect(df.col('formula')!.getTag('jmp.formula'), ':ints * 2');
    expect(df.col('charconstwidth')!.getTag(DG.TAGS.DESCRIPTION), 'Letters');
    expectArray(df.col('charconstwidth')!.getCategoryOrder(), ['d', 'c', 'b', 'a']);
    expect(JSON.parse(df.col('charconstwidth')!.getTag('jmp.properties')!).length, 2);
  });

  test('value labels show as the labels', async () => {
    const t = await fixture('example1');
    const ints = t.columns.find((c) => c.info.name === 'ints')!;
    const order = 'Value Order({Custom Order({3, 1})})';
    ints.properties = {
      properties: [{name: 'Value Order', text: order, node: parseJsl(order) as JslCall}],
      valueLabels: {keys: [1, 3], labels: ['One', 'Three']},
    };
    const col = jmpTableToDataFrame(t).col('ints')!;
    expect(col.type, DG.COLUMN_TYPE.STRING);
    expectArray(col.toList(), ['One', '2', 'Three', '4']);
    expectArray(col.getCategoryOrder(), ['Three', 'One']);
    expect(col.getTag('jmp.valueLabels'), '[[1,"One"],[3,"Three"]]');
  });

  test('every fixture converts', async () => {
    for (const name of FIXTURES) {
      const t = await fixture(name);
      const df = jmpTableToDataFrame(t, name);
      expect(df.rowCount, t.info.nrows, name);
    }
  });

  test('file handler', async () => {
    const bytes = await fixtureBytes('example1');
    const result: DG.DataFrame[] = await grok.functions.call('JmpReader:importJmp', {bytes});
    expect(result.length, 1);
    expect(result[0].rowCount, 4);
    expect(result[0].columns.length, 12);
  });

  // Through the platform, so the test bundle stays free of u2 (which the Node pass has no runtime for).
  test('preview builds and disposes', async () => {
    const u2 = (DG as any).U2;
    const before: number = u2.Scope.liveCount;
    const file = (await _package.files.list('tests')).find((f) => f.name === 'unsupported-column.jmp')!;
    const view: DG.ViewBase = await grok.functions.call('JmpReader:previewJmp', {file});
    expect(view.name, 'unsupported-column.jmp');
    expect(view.root.querySelector('[data-u2-name="summary"]') !== null, true);
    expectArray(Array.from(view.root.querySelectorAll('[data-u2-item]')).map((e) => (e as HTMLElement).dataset.u2Item),
      ['formula']);
    grok.shell.addView(view);
    view.close();
    await DG.delay(200);
    expect(u2.Scope.liveCount <= before, true, `live scopes ${u2.Scope.liveCount}, before ${before}`);
  });

  test('preview lists saved analyses', async () => {
    const file = (await _package.files.list('tests')).find((f) => f.name === 'saved-fit-model.jmp')!;
    const view: DG.ViewBase = await grok.functions.call('JmpReader:previewJmp', {file});
    try {
      const item = view.root.querySelector('[data-u2-item="Fit A"]') as HTMLElement | null;
      expect(item !== null, true);
      expect(/Not converted: JMP's Fit Model/.test(item!.textContent ?? ''), true, item!.textContent ?? '');
    } finally {
      view.close();
    }
  });
});

// The §8 target: a 100 MB file opens in under 10 s, and the main thread isn't blocked for more than 2 s at
// a time. The file is made in the test (large-jmp.ts), so it runs only with `grok test --benchmark`.
category('JMP: performance', () => {
  test('100 MB file', async () => {
    const bytes = await makeLargeJmp(await fixtureBytes('compressed'), rowsForSize(100e6));
    let last = performance.now();
    let longest = 0;
    const ticker = setInterval(() => {
      const now = performance.now();
      longest = Math.max(longest, now - last);
      last = now;
    }, 20);
    try {
      const t0 = performance.now();
      const table = await readJmp(bytes);
      const t1 = performance.now();
      const longestRead = longest;
      const df = jmpTableToDataFrame(table, 'large');
      const t2 = performance.now();
      const tv = grok.shell.addTableView(df);
      await DG.delay(100);
      const t3 = performance.now();
      tv.close();
      const s = (ms: number): string => (ms / 1000).toFixed(2);
      const report = `${(bytes.length / 1e6).toFixed(1)} MB, ${df.rowCount.toLocaleString()} rows: ` +
        `read ${s(t1 - t0)} s, DataFrame ${s(t2 - t1)} s, view ${s(t3 - t2)} s, total ${s(t3 - t0)} s; ` +
        `longest main-thread block ${s(longestRead)} s while reading, ${s(longest)} s overall`;
      console.log(`JMP performance: ${report}`);
      expect(t3 - t0 < 10000, true, report);
      return report;
    } finally {
      clearInterval(ticker);
    }
  }, {benchmark: true, timeout: 300000,
    skipReason: DG.Test.isInBenchmark ? undefined : 'run with grok test --benchmark'});
});
