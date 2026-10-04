import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

// The limits are the 5 s goal for applying a big layout with formula columns and 3-5x the quiet-machine times
// for the rest; smaller slowdowns show on the benchmark dashboard.
category('Benchmarks: Calculated columns and layouts', () => {
  const formulas = makeFormulas(10);
  let rows: number;
  let layout: string;
  let otherLayout: string;
  let df: DG.DataFrame | null = null;
  let secondView: DG.TableView | null = null;

  before(async () => {
    rows = DG.Test.isInBenchmark ? 50000 : 1000;
    const small = makeTable(20);
    for (const f of formulas)
      await small.columns.addNewCalculated(f.name, f.formula);
    const view = grok.shell.addTableView(small);
    layout = view.saveLayout().toJson();
    view.grid.columns.setOrder(small.columns.names().reverse());
    view.grid.props.rowHeight = 30;
    otherLayout = view.saveLayout().toJson();
    view.close();
    grok.shell.closeTable(small);
  });

  after(async () => {
    grok.shell.closeAll();
    df = null;
    secondView = null;
  });

  async function applyLayout(): Promise<number> {
    const table = makeTable(rows);
    const view = grok.shell.addTableView(table);
    const columnCount = table.columns.length + formulas.length;
    const start = performance.now();
    view.loadLayout(DG.ViewLayout.fromJson(layout));
    await awaitCheck(() => table.columns.length === columnCount, 'formula columns were not added', 300000, 20);
    df = table;
    return await waitIdle(start);
  }

  async function tableWithFormulas(): Promise<DG.DataFrame> {
    if (df == null)
      await applyLayout();
    return df!;
  }

  async function openSecondView(): Promise<number> {
    const table = await tableWithFormulas();
    const start = performance.now();
    secondView = grok.shell.addTableView(table);
    secondView.loadLayout(DG.ViewLayout.fromJson(layout));
    return await waitIdle(start);
  }

  test('Apply a layout with 200 formula columns', async () => expectFaster(await applyLayout(), 5000),
    {benchmark: true, timeout: 120000});

  test('Open a second view with the layout', async () => expectFaster(await openSecondView(), 500), {benchmark: true, timeout: 120000});

  test('Switch a view to another layout', async () => {
    if (secondView == null)
      await openSecondView();
    const start = performance.now();
    secondView!.loadLayout(DG.ViewLayout.fromJson(otherLayout));
    return expectFaster(await waitIdle(start), 300);
  }, {benchmark: true, timeout: 120000});

  test('Add a formula column', async () => {
    const table = await tableWithFormulas();
    const start = performance.now();
    const column = await table.columns.addNewCalculated('added', '${a} * 2 + ${f9_1}');
    const ms = performance.now() - start;
    table.columns.remove(column.name);
    return expectFaster(ms, 200);
  }, {benchmark: true, timeout: 120000});

  test('Edit a cell that 10 levels of formulas depend on', async () => {
    const table = await tableWithFormulas();
    const deepest = table.col('f9_13')!;
    let row = 0;
    while (row < table.rowCount && deepest.isNone(row))
      row++;
    if (row === table.rowCount)
      throw new Error('f9_13 has no values');
    const before = deepest.get(row);
    const start = performance.now();
    table.set('b', row, table.get('b', row) + 10);
    await awaitCheck(() => deepest.get(row) !== before, 'dependent formulas were not recalculated', 60000, 5);
    return expectFaster(performance.now() - start, 500);
  }, {benchmark: true, timeout: 120000});

  test('Calculate 200 formulas', async () => {
    const table = makeTable(rows);
    const start = performance.now();
    for (const f of formulas)
      await table.columns.addNewCalculated(f.name, f.formula);
    return expectFaster(performance.now() - start, 5000);
  }, {benchmark: true, timeout: 120000});

  test('Formulas failing on every row', async () => {
    const table = makeTable(rows);
    const start = performance.now();
    for (const formula of ['RegExpReplace(${s}, "(", "x")', 'ParseInt(${s})']) {
      const column = await table.columns.addNewCalculated('failing', formula);
      expect(column.stats.missingValueCount, rows);
      table.columns.remove(column.name);
    }
    return expectFaster(performance.now() - start, 2500);
  }, {benchmark: true, timeout: 120000});
}, {owner: 'dkovalyov@datagrok.ai', clear: false});

function expectFaster(ms: number, benchmarkLimitMs: number): string {
  if (DG.Test.isInBenchmark && ms > benchmarkLimitMs)
    throw new Error(`${Math.round(ms)} ms, expected under ${benchmarkLimitMs} ms`);
  return `${Math.round(ms)} ms`;
}

/** Resolves once the page has had no task longer than 50 ms for [quietMs], or after [maxMs]; returns when it was last busy. */
async function waitIdle(start: number, quietMs: number = 200, maxMs: number = 60000): Promise<number> {
  let lastBusy = performance.now();
  while (performance.now() - lastBusy < quietMs && performance.now() - start < maxMs) {
    const t = performance.now();
    await delay(20);
    if (performance.now() - t > 70)
      lastBusy = performance.now();
  }
  return Math.round(lastBusy - start);
}

function makeTable(rowCount: number): DG.DataFrame {
  let seed = 42;
  const random = () => (seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff) / 0x7fffffff;
  const floats = (name: string, emptyShare: number) => DG.Column.fromFloat32Array(name, Float32Array.from(
    {length: rowCount}, () => random() < emptyShare ? DG.FLOAT_NULL : Math.round(random() * 10000) / 100));
  const ints = (name: string, emptyShare: number) => DG.Column.fromInt32Array(name,
    Int32Array.from({length: rowCount}, () => random() < emptyShare ? DG.INT_NULL : (random() * 1000) | 0));
  const categories = (name: string, values: string[]) => DG.Column.fromIndexes(name, values,
    Int32Array.from({length: rowCount}, () => (random() * values.length) | 0));
  const flags = Array.from({length: rowCount}, () => random() < 0.5);
  const dates = Array.from({length: rowCount},
    () => random() < 0.75 ? null : dayjs(Date.UTC(2015, 0, 1) + random() * 3e11));
  return DG.DataFrame.fromColumns([
    floats('a', 0.1), floats('b', 0.1), floats('c', 0.3), ints('n1', 0.2), ints('n2', 0.2),
    categories('cat', ['', 'c1', 'c2', 'c3', 'c4', 'c5']),
    categories('s', ['', 't.e)st.a', 'x.y.z', 'alpha.beta', 'HT-Sol (2)']),
    DG.Column.fromList(DG.TYPE.BOOL, 'flag', flags), DG.Column.fromList(DG.TYPE.DATE_TIME, 'd', dates)]);
}

/** [rounds] x 20 formulas mixing the constructs of real layouts, chained across rounds, plus one that always fails. */
function makeFormulas(rounds: number): {name: string, formula: string}[] {
  const col = (name: string) => '${' + name + '}';
  const [a, b, c, n1, n2, cat, s, flag, d] = ['a', 'b', 'c', 'n1', 'n2', 'cat', 's', 'flag', 'd'].map(col);
  const list: {name: string, formula: string}[] = [];
  for (let r = 0; r < rounds; r++) {
    const k = r + 2;
    const prev = (i: number) => col(`f${r - 1}_${i}`);
    [
      `${a} * ${k} + ${b}`,
      `(${a} - ${c}) / (${b} + ${k})`,
      `if(${cat} == "c1" || ${cat} == "c2", ${a} * ${k}, if(${cat} == "c3", ${b}, null))`,
      `if(${n1} != null, ${n1} * 2, ${n2})`,
      `${cat} + " (" + ToString(Round10(${a}, 2)) + ")"`,
      `ToUpperCase(${s})`,
      `SplitString(${s}, ".", 1)`,
      `Year(${d})`,
      `InDays(DateDiff(${d}, Date(2020, 1, ${k})))`,
      `Qnum(${a} * ${k}, if(${a} > 50, ">", "="))`,
      `if(Qualifier(${col(`f${r}_9`)}) == ">", "high", "normal")`,
      `BinBySpecificLimits(${a}, [10, 50, 90])`,
      `RegExpReplace(${s}, "[.]", "-")`,
      r === 0 ? `Round10(${b} / (Abs(${c}) + 1), 3)` : `Round10(${prev(13)} * 1.01 + ${prev(0)}, 3)`,
      r === 0 ? `${a} + ${n1}` : `${prev(1)} + ${prev(3)}`,
      `case when ${a} > 75 then "A" when ${a} > 50 then "B" when ${a} > 25 then "C" else "D" end`,
      `${flag} && ${a} > ${k * 5}`,
      `Max([${a}, ${b}, ${c}])`,
      `${cat} in ["c1", "c3", "c5"]`,
      `if(IsEmpty(${d}), "no date", "dated")`,
    ].forEach((formula, i) => list.push({name: `f${r}_${i}`, formula}));
  }
  list.push({name: 'failing', formula: `RegExpReplace(${s}, "(", "x")`});
  return list;
}
