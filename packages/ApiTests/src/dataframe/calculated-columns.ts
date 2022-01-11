import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { Subscription } from 'rxjs';
import { after, category, expect, test } from '../test';


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
      let column = await df.columns.addNewCalculated('new', '${x}+${y}-${z}');
      expect(df.columns.contains(column.name), true);
      expect(column.tags[DG.TAGS.FORMULA], '${x}+${y}-${z}');
      expect(column.get(0), -2);
      expect(column.get(1), -1);
      expect(column.get(2), 0);
    } finally {
      df.columns.remove('new');
    }
  });

  test('Add new column dialog', async () => {
    let tv: DG.TableView;
    subs.push(grok.events.onDialogShown.subscribe((d) => {
      expect(d.title, 'Add New Column');
      dialogs.push(d);
    }));
    try {
      tv = grok.shell.addTableView(df);
      df.dialogs.addNewColumn();
    } finally {
      tv!.close();
      grok.shell.closeTable(df);
    }
  });

  test('Edit formula dialog', async () => {
    subs.push(grok.events.onDialogShown.subscribe((d) => {
      expect(d.title, 'Edit Column Formula');
      dialogs.push(d);
    }));
    try {
      let column = await df.columns.addNewCalculated('editable', '0');
      column.dialogs.editFormula();
    } finally {
      df.columns.remove('editable');
    }
  });

  test('Calculated columns addition event', async () => {
    const t = df.clone();
    subs.push(t.onColumnsAdded.subscribe((data) => data.args.columns.forEach((column: DG.Column) => {
      expect(column.name, column.tags.has(DG.TAGS.FORMULA) ? 'calculated column' : 'regular column');
    })));
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    t.columns.addNewInt('regular column').init(1);
  });

  test('Calculated columns deletion event', async () => {
    const t = df.clone();
    subs.push(t.onColumnsRemoved.subscribe((data) => data.args.columns.forEach((column: DG.Column) => {
      expect(column.name, column.tags.has(DG.TAGS.FORMULA) ? 'calculated column' : 'regular column');
    })));
    await t.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
    t.columns.addNewInt('regular column').init(1);
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    dialogs.forEach((d) => d.close());
  });
});
