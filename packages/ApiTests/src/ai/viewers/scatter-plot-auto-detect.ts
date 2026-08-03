import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, look} from '../helpers';

// ScatterPlot X/Y + whisker auto-detection (ScatterPlotLook.auto).

function nums(name: string, vals: number[]): [string, string, number[]] {
  return [name, 'double', vals];
}

function expectAxisPair(c: DG.ScatterPlotViewer, candidates: [string, string]): {xBase: string, yBase: string} {
  const xBase = c.props['xColumnName'] as string;
  const yBase = c.props['yColumnName'] as string;
  expect(xBase === candidates[0] || xBase === candidates[1], true);
  expect(yBase === candidates[0] || yBase === candidates[1], true);
  expect(xBase !== yBase, true);
  return {xBase, yBase};
}

function expectCleared(value: any): void {
  expect(value == null || value === '', true);
}

category('AI: Viewers: ScatterPlot auto-detect', () => {

  // === X/Y auto-pick (no whiskers) ===

  test('two numerical columns pick distinct X and Y', async () => {
    const c = df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]), nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    expectAxisPair(c, ['alpha', 'beta']);
  });

  test('single numerical column places same column on both axes', async () => {
    const c = df([nums('only', [1, 2, 3, 4, 5, 6, 7, 8])]).plot.scatter() as DG.ScatterPlotViewer;
    expect(c.props['xColumnName'], 'only');
    expect(c.props['yColumnName'], 'only');
  });

  test('numerical preferred over categorical when both present', async () => {
    const c = df([
      ['cat', 'string', ['a', 'b', 'c', 'a', 'b', 'c', 'a', 'b']],
      nums('num', [1, 2, 3, 4, 5, 6, 7, 8]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    expect(c.props['xColumnName'], 'num');
    expect(c.props['yColumnName'], 'cat');
  });

  test('categorical-only dataframe falls back to categoricals', async () => {
    const c = df([
      ['cat1', 'string', ['a', 'b', 'c', 'a', 'b', 'c']],
      ['cat2', 'string', ['x', 'y', 'z', 'x', 'y', 'z']],
    ]).plot.scatter() as DG.ScatterPlotViewer;
    expectAxisPair(c, ['cat1', 'cat2']);
  });

  test('tilde-prefixed columns excluded from auto-pick when others exist', async () => {
    const c = df([
      nums('~hidden1', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('~hidden2', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('visible1', [100, 200, 300, 400, 500, 600, 700, 800]),
      nums('visible2', [9, 18, 27, 36, 45, 54, 63, 72]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    const {xBase, yBase} = expectAxisPair(c, ['visible1', 'visible2']);
    expect(xBase.startsWith('~'), false);
    expect(yBase.startsWith('~'), false);
  });

  // === Whisker auto-detect: suffix family coverage ===

  function buildPair(suffix: string): DG.DataFrame {
    return df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a' + suffix, [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('b', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('b' + suffix, [9, 18, 28, 38, 48, 58, 68, 78]),
    ]);
  }

  function buildMinMax(suffixMin: string, suffixMax: string): DG.DataFrame {
    return df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a' + suffixMin, [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('a' + suffixMax, [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      nums('b', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('b' + suffixMin, [9, 18, 28, 38, 48, 58, 68, 78]),
      nums('b' + suffixMax, [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
  }

  for (const s of [' whisker', ' error']) {
    const suffix = s;
    test(`suffix "${suffix} min" / "${suffix} max" auto-detected`, async () => {
      const c = buildMinMax(suffix + ' min', suffix + ' max').plot.scatter() as DG.ScatterPlotViewer;
      const {xBase, yBase} = expectAxisPair(c, ['a', 'b']);
      expect(c.props['xWhiskerMinColumnName'], xBase + suffix + ' min');
      expect(c.props['xWhiskerMaxColumnName'], xBase + suffix + ' max');
      expect(c.props['yWhiskerMinColumnName'], yBase + suffix + ' min');
      expect(c.props['yWhiskerMaxColumnName'], yBase + suffix + ' max');
    });
  }

  for (const r of [' whisker range', ' error range']) {
    const suffix = r;
    test(`suffix "${suffix}" auto-detected (no min/max)`, async () => {
      const c = buildPair(suffix).plot.scatter() as DG.ScatterPlotViewer;
      const {xBase, yBase} = expectAxisPair(c, ['a', 'b']);
      expect(c.props['xWhiskerRangeColumnName'], xBase + suffix);
      expect(c.props['yWhiskerRangeColumnName'], yBase + suffix);
      expectCleared(c.props['xWhiskerMinColumnName']);
      expectCleared(c.props['yWhiskerMinColumnName']);
    });
  }

  for (const sd of [' std', ' stddev']) {
    const suffix = sd;
    test(`bare "${suffix}" suffix folded into range`, async () => {
      const c = buildPair(suffix).plot.scatter() as DG.ScatterPlotViewer;
      const {xBase, yBase} = expectAxisPair(c, ['a', 'b']);
      expect(c.props['xWhiskerRangeColumnName'], xBase + suffix);
      expect(c.props['yWhiskerRangeColumnName'], yBase + suffix);
      expectCleared(c.props['xWhiskerMinColumnName']);
      expectCleared(c.props['xWhiskerMaxColumnName']);
    });
  }

  test('range coexists with min/max when all three suffixes present', async () => {
    const c = df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a error range', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      nums('a error min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('a error max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      nums('b', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('b error range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
      nums('b error min', [9, 18, 28, 38, 48, 58, 68, 78]),
      nums('b error max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    const {xBase, yBase} = expectAxisPair(c, ['a', 'b']);
    expect(c.props['xWhiskerRangeColumnName'], xBase + ' error range');
    expect(c.props['xWhiskerMinColumnName'], xBase + ' error min');
    expect(c.props['xWhiskerMaxColumnName'], xBase + ' error max');
    expect(c.props['yWhiskerRangeColumnName'], yBase + ' error range');
    expect(c.props['yWhiskerMinColumnName'], yBase + ' error min');
    expect(c.props['yWhiskerMaxColumnName'], yBase + ' error max');
  });

  test('min-only suffix auto-detected without a matching max', async () => {
    const c = buildPair(' whisker min').plot.scatter() as DG.ScatterPlotViewer;
    const {xBase, yBase} = expectAxisPair(c, ['a', 'b']);
    expect(c.props['xWhiskerMinColumnName'], xBase + ' whisker min');
    expect(c.props['yWhiskerMinColumnName'], yBase + ' whisker min');
    expectCleared(c.props['xWhiskerMaxColumnName']);
    expectCleared(c.props['xWhiskerRangeColumnName']);
  });

  test('fallback axis pick skips whisker columns', async () => {
    const c = df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('a max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      nums('plain', [100, 200, 300, 400, 500, 600, 700, 800]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    expect(c.props['xColumnName'], 'a');
    expect(c.props['xWhiskerMinColumnName'], 'a min');
    expect(c.props['xWhiskerMaxColumnName'], 'a max');
    expect(c.props['yColumnName'], 'plain');
  });

  // === Tag-based detection ===

  test('tag "whisker.range" on numeric column drives axis pick + range slot', async () => {
    const d = df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('rng1', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      nums('rng2', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    d.col('alpha')!.tags['whisker.range'] = 'rng1';
    d.col('beta')!.tags['whisker.range'] = 'rng2';
    const c = d.plot.scatter() as DG.ScatterPlotViewer;
    const {xBase, yBase} = expectAxisPair(c, ['alpha', 'beta']);
    const expRange = (axis: string) => axis === 'alpha' ? 'rng1' : 'rng2';
    expect(c.props['xWhiskerRangeColumnName'], expRange(xBase));
    expect(c.props['yWhiskerRangeColumnName'], expRange(yBase));
  });

  // === Per-axis re-detect after options apply ===
  // Passing {x: <base>} re-runs the per-axis whisker detector, so whiskers still light up.

  function expectAxisWhiskers(axis: 'x' | 'y', base: string): void {
    const opts: {[k: string]: string} = {};
    opts[axis] = base;
    const c = buildMinMax(' min', ' max').plot.scatter(opts) as DG.ScatterPlotViewer;
    const p = c.props as any;
    expect(p[axis + 'ColumnName'], base);
    expect(p[axis + 'WhiskerMinColumnName'], base + ' min');
    expect(p[axis + 'WhiskerMaxColumnName'], base + ' max');
  }

  test('per-axis whisker re-detect honors xColumnName passed via options', async () =>
    expectAxisWhiskers('x', 'a'));

  test('per-axis whisker re-detect honors yColumnName passed via options', async () =>
    expectAxisWhiskers('y', 'b'));

}, {owner: 'agolovko@datagrok.ai'});

// Same auto-detection, exercised through a trellis inner viewer.

const splitRow: [string, string, string[]] = ['split', 'string', ['s', 't', 's', 't', 's', 't', 's', 't']];

function trellisInner(d: DG.DataFrame): {[k: string]: any} {
  const tv = DG.Viewer.trellisPlot(d, {viewerType: DG.VIEWER.SCATTER_PLOT, xColumnNames: ['split']});
  return look(tv)['innerViewerLook'];
}

function expectInnerAxisPair(inner: {[k: string]: any}, candidates: [string, string]): {xBase: string, yBase: string} {
  const xBase = inner['xColumnName'];
  const yBase = inner['yColumnName'];
  expect(xBase === candidates[0] || xBase === candidates[1], true);
  expect(yBase === candidates[0] || yBase === candidates[1], true);
  expect(xBase !== yBase, true);
  return {xBase, yBase};
}

category('AI: Viewers: ScatterPlot auto-detect (as trellis inner viewer)', () => {

  test('two numerical columns pick distinct X and Y (trellis inner)', async () => {
    const inner = trellisInner(df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]), nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]), splitRow,
    ]));
    expectInnerAxisPair(inner, ['alpha', 'beta']);
  });

  test('tilde-prefixed columns excluded from auto-pick (trellis inner)', async () => {
    const inner = trellisInner(df([
      nums('~hidden1', [1, 2, 3, 4, 5, 6, 7, 8]), nums('~hidden2', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('visible1', [100, 200, 300, 400, 500, 600, 700, 800]), nums('visible2', [9, 18, 27, 36, 45, 54, 63, 72]),
      splitRow,
    ]));
    const {xBase, yBase} = expectInnerAxisPair(inner, ['visible1', 'visible2']);
    expect(xBase.startsWith('~'), false);
    expect(yBase.startsWith('~'), false);
  });

  test('suffix " whisker min" / " whisker max" auto-detected (trellis inner)', async () => {
    const inner = trellisInner(df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a whisker min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('a whisker max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      nums('b', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('b whisker min', [9, 18, 28, 38, 48, 58, 68, 78]),
      nums('b whisker max', [11, 22, 32, 42, 52, 62, 72, 82]),
      splitRow,
    ]));
    const {xBase, yBase} = expectInnerAxisPair(inner, ['a', 'b']);
    expect(inner['xWhiskerMinColumnName'], xBase + ' whisker min');
    expect(inner['xWhiskerMaxColumnName'], xBase + ' whisker max');
    expect(inner['yWhiskerMinColumnName'], yBase + ' whisker min');
    expect(inner['yWhiskerMaxColumnName'], yBase + ' whisker max');
  });

  test('suffix " error range" auto-detected (no min/max, trellis inner)', async () => {
    const inner = trellisInner(df([
      nums('a', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('a error range', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      nums('b', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('b error range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
      splitRow,
    ]));
    const {xBase, yBase} = expectInnerAxisPair(inner, ['a', 'b']);
    expect(inner['xWhiskerRangeColumnName'], xBase + ' error range');
    expect(inner['yWhiskerRangeColumnName'], yBase + ' error range');
  });

  test('tag "whisker.range" drives axis pick + range slot (trellis inner)', async () => {
    const d = df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('rng1', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      nums('rng2', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
      splitRow,
    ]);
    d.col('alpha')!.tags['whisker.range'] = 'rng1';
    d.col('beta')!.tags['whisker.range'] = 'rng2';
    const inner = trellisInner(d);
    const {xBase, yBase} = expectInnerAxisPair(inner, ['alpha', 'beta']);
    const expRange = (axis: string) => axis === 'alpha' ? 'rng1' : 'rng2';
    expect(inner['xWhiskerRangeColumnName'], expRange(xBase));
    expect(inner['yWhiskerRangeColumnName'], expRange(yBase));
  });

}, {owner: 'agolovko@datagrok.ai'});
