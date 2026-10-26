import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import dayjs from 'dayjs';
import {after, category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';


category('DataFrame: Calculated columns', () => {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9]),
  ]);
  const subs: Subscription[] = [];
  const dialogs: DG.Dialog[] = [];

  test('Create a calculated column', async () => {
    try {
      const column = await df.columns.addNewCalculated('new', '${x}+${y}-${z}');
      expect(df.columns.contains(column.name), true);
      expect(column.meta.formula, '${x}+${y}-${z}');
      expect(column.get(0), -2);
      expect(column.get(1), -1);
      expect(column.get(2), 0);
    } finally {
      df.columns.remove('new');
    }
  });

  test('Create a calculated column with script formula', async () => {
    try {
      const column = await df.columns.addNewCalculated('new', 'ApiTests:FormulaScript(${x})');
      expect(df.columns.contains(column.name), true);
      expect(column.get(0), 11);
    } finally {
      df.columns.remove('new');
    }
  }, {skipReason: typeof process !== 'undefined' ? 'client package functions are not loaded in NodeJS' : undefined});

  test('Create a calculated column with async formula', async () => {
    try {
      const column = await df.columns.addNewCalculated('new', 'ApiTests:testIntAsync(${x})');
      expect(df.columns.contains(column.name), true);
      expect(column.get(0), 11);
    } finally {
      df.columns.remove('new');
    }
  }, {skipReason: typeof process !== 'undefined' ? 'client package functions are not loaded in NodeJS' : undefined});

  test('Rows the formula fails on', async () => {
    const t = DG.DataFrame.fromColumns([DG.Column.fromStrings('s', ['2023-01-05', 'n/a', '2023-03-07'])]);
    const empty = await t.columns.addNewCalculated('empty', 'DateParse(${s})');
    expect(empty.isNone(1), true);
    expect(empty.getTag(DG.Tags.FormulaErrorBehavior) == null, true);

    const filled = await t.columns.addNewCalculated('filled', 'DateParse(${s})',
      {type: 'datetime', onError: {mode: 'value', value: dayjs.utc('1900-01-01'), errorColumn: true}});
    expect(filled.get(1)!.valueOf(), dayjs.utc('1900-01-01').valueOf());
    const errors = t.col('filled errors')!;
    expect(errors.get(1).includes('n/a'), true);
    expect(errors.isNone(0), true);
    expect(JSON.parse(filled.getTag(DG.Tags.FormulaErrorBehavior)!).errorColName, 'filled errors');

    let rejected = false;
    try {
      await t.columns.addNewCalculated('strict', 'DateParse(${s})', {onError: {mode: 'stop'}});
    } catch (_) {
      rejected = true;
    }
    expect(rejected, true);
    expect(t.col('strict'), null);

    const q = DG.DataFrame.fromColumns([DG.Column.fromStrings('s', ['1.5', 'n/a'])]);
    const qnum = await q.columns.addNewCalculated('q', 'ParseFloat(${s})',
      {type: 'qnum', onError: {mode: 'value', value: DG.Qnum.less(5)}});
    expect(DG.Qnum.qualifier(qnum.get(1)), '<');
    expect(DG.Qnum.getValue(qnum.get(1)), 5);
  });

  test('Add new column dialog', async () => {
    if ((await grok.dapi.packages.filter('PowerPack').list({pageSize: 5})).length > 0)
      return 'Skipped because PowerPack is installed';
    let shown = false;
    subs.push(grok.events.onDialogShown.subscribe((d: any) => {
      if (d.title == 'Add New Column')
        shown = true;
      dialogs.push(d);
    }));
    const tv = grok.shell.addTableView(df);
    try {
      await df.dialogs.addNewColumn();
      await waitUntil(() => shown, 'Dialog not found');
    } finally {
      tv.close();
      grok.shell.closeTable(df);
    }
  });

  test('Edit formula dialog', async () => {
    if ((await grok.dapi.packages.filter('PowerPack').list({pageSize: 5})).length > 0)
      return 'Skipped because PowerPack is installed';
    let shown = false;
    subs.push(grok.events.onDialogShown.subscribe((d: DG.Dialog) => {
      if (d.title == 'Add New Column')
        shown = true;
      dialogs.push(d);
    }));
    try {
      const column = await df.columns.addNewCalculated('editable', '0');
      column.meta.dialogs.editFormula();
      await waitUntil(() => shown, 'Dialog not found');
    } finally {
      df.columns.remove('editable');
    }
  });

  test('Calculated columns addition event', async () => {
    const t = df.clone();
    let added = false;
    subs.push(t.onColumnsAdded.subscribe((data: any) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.meta.formula !== null && column.name === 'calculated column')
          added = true;
      })));

    t.columns.addNewInt('regular column').init(1);
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    await waitUntil(() => added, 'Failed to add a calculated column');
  });

  test('Calculated columns deletion event', async () => {
    const t = df.clone();
    let removed = false;
    subs.push(t.onColumnsRemoved.subscribe((data: any) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.meta.formula !== null && column.name === 'calculated column')
          removed = true;
      })));

    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    t.columns.remove('calculated column');
    await waitUntil(() => removed, 'Failed to delete a calculated column');
  });

  test('A formula naming its column as a string recalculates after an edit', async () => {
    const t = DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3])]);
    await t.columns.addNewCalculated('ctrl', 'Avg($[x])');
    await t.columns.addNewCalculated('byName', 'Avg(Column("x"))');
    await t.columns.addNewCalculated('first', 'Value("x", 0)');
    t.col('x')!.set(0, 10);
    await waitUntil(() => t.col('ctrl')!.get(0) === 5, 'Avg($[x]) did not recalculate');
    expect(t.col('byName')!.get(0), 5);
    expect(t.col('first')!.get(0), 10);
  }, {skipReason: 'core#1128'});

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    dialogs.forEach((d) => d.close());
  });
}, {owner: 'mdolotova@datagrok.ai'});
