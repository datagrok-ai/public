import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {addSavedAnalyses, jmpTableToDataFrame, savedAnalyses, skippedColumns} from '../formats/jmp/jmp-to-dataframe';
import {readJmp} from '../formats/jmp/jmp-parser';
import {fixture, fixtureBytes, FIXTURES} from './jmp-reader-tests';

category('JMP: dataframe', () => {
  test('types and nulls', async () => {
    const df = jmpTableToDataFrame(await fixture('example1'), 'example1');
    expect(df.rowCount, 4);
    expect(df.col('ints')!.type, DG.COLUMN_TYPE.INT);
    expect(df.col('floats')!.type, DG.COLUMN_TYPE.FLOAT);
    expect(df.col('charconstwidth')!.type, DG.COLUMN_TYPE.STRING);
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

  test('every fixture converts', async () => {
    for (const name of FIXTURES) {
      const t = await fixture(name);
      const df = jmpTableToDataFrame(t, name);
      expect(df.rowCount, t.info.nrows, name);
    }
  });

  test('file handler', async () => {
    const result: DG.DataFrame[] = await grok.functions.call('Jmp:importJmp', {bytes: await fixtureBytes('example1')});
    expect(result.length, 1);
    expect(result[0].rowCount, 4);
    expect(result[0].columns.length, 12);
  });

  // Through the platform, so the test bundle stays free of u2 (which the Node pass has no runtime for).
  test('preview builds and disposes', async () => {
    const u2 = (DG as any).U2;
    const before: number = u2.Scope.liveCount;
    const file = (await _package.files.list('tests')).find((f) => f.name === 'unsupported-column.jmp')!;
    const view: DG.ViewBase = await grok.functions.call('Jmp:previewJmp', {file});
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
    const view: DG.ViewBase = await grok.functions.call('Jmp:previewJmp', {file});
    try {
      const item = view.root.querySelector('[data-u2-item="Fit A"]') as HTMLElement | null;
      expect(item !== null, true);
      expect(/Not converted: JMP's Fit Model/.test(item!.textContent ?? ''), true, item!.textContent ?? '');
    } finally {
      view.close();
    }
  });
});
