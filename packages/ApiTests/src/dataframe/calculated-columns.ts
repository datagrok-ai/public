import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';


category('DataFrame: Calculated columns', () => {
  async function testActionJoinSelectiveColumns(funcName: string): Promise<void> {
    const testDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.INT, 'col1', [1, 2, 3]),
      DG.Column.fromList(DG.TYPE.STRING, 'col2', ['a', 'b', 'c']),
      DG.Column.fromList(DG.TYPE.FLOAT, 'col3', [1.5, 2.5, 3.5]),
    ]);
    const func = DG.Func.find({package: 'DevTools', name: funcName})[0];
    if (!func)
      throw new Error(`${funcName} not found in DevTools package`);
    const call = func.prepare({
      data: testDf,
      col1: testDf.col('col1'),
      col2: testDf.col('col2'),
      col3: testDf.col('col3'),
      out: ['joinedCol1', 'joinedCol3'],
    });
    await call.call();

    expect(testDf.columns.contains('joinedCol1'), true, 'joinedCol1 should be added');
    expect(testDf.columns.contains('joinedCol2'), false, 'joinedCol2 should NOT be added');
    expect(testDf.columns.contains('joinedCol3'), true, 'joinedCol3 should be added');
  }

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
  });

  test('Create a calculated column with async formula', async () => {
    try {
      const column = await df.columns.addNewCalculated('new', 'ApiTests:testIntAsync(${x})');
      expect(df.columns.contains(column.name), true);
      expect(column.get(0), 11);
    } finally {
      df.columns.remove('new');
    }
  });

  test('Add new column dialog', () => new Promise(async (resolve, reject) => {
    if ((await grok.dapi.packages.filter('PowerPack').list({pageSize: 5})).length > 0)
      resolve('Skipped because PowerPack is installed');
    else {
      let tv: DG.TableView;
      subs.push(grok.events.onDialogShown.subscribe((d: any) => {
        if (d.title == 'Add New Column')
          resolve('OK');
        dialogs.push(d);
      }));
      setTimeout(() => {
        // eslint-disable-next-line prefer-promise-reject-errors
        reject('Dialog not found');
      }, 1000);
      try {
        tv = grok.shell.addTableView(df);
        await df.dialogs.addNewColumn();
      } finally {
        tv!.close();
        grok.shell.closeTable(df);
      }
    }
  }));

  test('Edit formula dialog', () => new Promise(async (resolve, reject) => {
    if ((await grok.dapi.packages.filter('PowerPack').list({pageSize: 5})).length > 0)
      resolve('Skipped because PowerPack is installed');
    else {
      subs.push(grok.events.onDialogShown.subscribe((d: DG.Dialog) => {
        if (d.title == 'Add New Column')
          resolve('OK');
        dialogs.push(d);
      }));
      try {
        setTimeout(() => {
          // eslint-disable-next-line prefer-promise-reject-errors
          reject('Dialog not found');
        }, 1000);
        const column = await df.columns.addNewCalculated('editable', '0');
        column.meta.dialogs.editFormula();
      } finally {
        df.columns.remove('editable');
      }
    }
  }));

  test('Calculated columns addition event', () => new Promise(async (resolve, reject) => {
    const t = df.clone();
    subs.push(t.onColumnsAdded.subscribe((data: any) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.meta.formula !== null && column.name === 'calculated column')
          resolve('OK');
      })));

    setTimeout(() => reject(new Error('Failed to add a calculated column')), 50);
    t.columns.addNewInt('regular column').init(1);
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
  }));

  test('Calculated columns deletion event', () => new Promise(async (resolve, reject) => {
    const t = df.clone();
    subs.push(t.onColumnsRemoved.subscribe((data: any) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.meta.formula !== null && column.name === 'calculated column')
          resolve('OK');
      })));

    setTimeout(() => reject(new Error('Failed to delete a calculated column')), 100);
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    t.columns.addNewInt('regular column').init(1);
    await DG.delay(50);
    t.columns.remove('regular column');
    t.columns.remove('calculated column');
  }));

  test('Action join function selective columns', async () => {
    await testActionJoinSelectiveColumns('testFunctionJoin');
  });

  test('Action join script selective columns', async () => {
    await testActionJoinSelectiveColumns('TestFunctionScriptJoin');
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    dialogs.forEach((d) => d.close());
  });
}, {owner: 'mdolotova@datagrok.ai'});
