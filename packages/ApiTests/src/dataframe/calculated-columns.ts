import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {after, category, expect, test} from '@datagrok-libraries/utils/src/test';


category('DataFrame', () => {
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
      expect(column.tags[DG.TAGS.FORMULA], '${x}+${y}-${z}');
      expect(column.get(0), -2);
      expect(column.get(1), -1);
      expect(column.get(2), 0);
    } finally {
      df.columns.remove('new');
    }
  });

  test('Add new column dialog', () => new Promise(async (resolve, reject) => {
    if ((await grok.dapi.packages.filter('PowerPack').list()).length > 0)
      resolve('Skipped because PowerPack is installed');
    else {
      let tv: DG.TableView;
      subs.push(grok.events.onDialogShown.subscribe((d) => {
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
    if ((await grok.dapi.packages.filter('PowerPack').list()).length > 0)
      resolve('Skipped because PowerPack is installed');
    else {
      subs.push(grok.events.onDialogShown.subscribe((d) => {
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
        column.dialogs.editFormula();
      } finally {
        df.columns.remove('editable');
      }
    }
  }));

  test('Calculated columns addition event', () => new Promise(async (resolve, reject) => {
    const t = df.clone();
    subs.push(t.onColumnsAdded.subscribe((data) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.tags.has(DG.TAGS.FORMULA) && column.name === 'calculated column')
          resolve('OK');
    })));

    setTimeout(() => reject('Failed to add a calculated column'), 50);
    t.columns.addNewInt('regular column').init(1);
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
  }));

  test('Calculated columns deletion event', () => new Promise(async (resolve, reject) => {
    const t = df.clone();
    subs.push(t.onColumnsRemoved.subscribe((data) =>
      data.args.columns.forEach((column: DG.Column) => {
        if (column.tags.has(DG.TAGS.FORMULA) && column.name === 'calculated column')
          resolve('OK');
    })));

    setTimeout(() => reject('Failed to delete a calculated column'), 50);
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    t.columns.addNewInt('regular column').init(1);
    t.columns.remove('regular column');
    t.columns.remove('calculated column');
  }));

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    dialogs.forEach((d) => d.close());
  });
});
