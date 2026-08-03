import { after, before, category, expect, test } from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Undo', () => {
  const nodeSkip = typeof process !== 'undefined' ? 'NodeJS environment' : undefined;
  let view: DG.TableView;
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(100);
    view = grok.shell.addTableView(df);
    DG.UndoService.clear();
  });

  after(async () => {
    DG.UndoService.clear();
    view?.close();
    grok.shell.closeTable(df);
  });

  test('push and undo', async () => {
    const column = df.col('age')!;
    const before = column.get(0);
    column.set(0, 100);
    DG.UndoService.push('Set age', () => column.set(0, before), {
      redo: () => column.set(0, 100),
      context: df,
    });

    expect(grok.shell.canUndo, true);
    expect(DG.UndoService.undoName, 'Set age');

    grok.shell.undo();
    expect(column.get(0), before);
    expect(grok.shell.canRedo, true);

    grok.shell.redo();
    expect(column.get(0), 100);

    grok.shell.undo();
    expect(column.get(0), before);
  }, { skipReason: nodeSkip, node: false });

  test('multi-level LIFO', async () => {
    DG.UndoService.clear();
    const column = df.col('age')!;
    const original = column.get(1);

    for (const value of [10, 20, 30]) {
      const previous = column.get(1);
      column.set(1, value);
      DG.UndoService.push(`Set ${value}`, () => column.set(1, previous), { context: df });
    }

    grok.shell.undo();
    expect(column.get(1), 20);
    grok.shell.undo();
    expect(column.get(1), 10);
    grok.shell.undo();
    expect(column.get(1), original);
    expect(grok.shell.canUndo, false);
  }, { skipReason: nodeSkip, node: false });

  test('onUndo fires', async () => {
    DG.UndoService.clear();
    let fired: string | null = null;
    const sub = grok.shell.onUndo.subscribe((args) => fired = args.name);
    try {
      DG.UndoService.push('Marker', () => {}, { context: df });
      grok.shell.undo();
      expect(fired, 'Marker');
    } finally {
      sub.unsubscribe();
    }
  }, { skipReason: nodeSkip, node: false });

  test('records are dropped when the table is closed', async () => {
    DG.UndoService.clear();
    const temp = grok.data.demo.demog(10);
    const tempView = grok.shell.addTableView(temp);
    DG.UndoService.push('Temp op', () => {}, { context: temp });
    expect(DG.UndoService.undoName, 'Temp op');

    tempView.close();
    grok.shell.closeTable(temp);
    // closing the view pushes its own record, so assert the table's record is gone by name
    expect(DG.UndoService.undoName !== 'Temp op', true);
  }, { skipReason: nodeSkip, node: false });

  test('no redo without a forward action', async () => {
    DG.UndoService.clear();
    DG.UndoService.push('One way', () => {}, { context: df });
    grok.shell.undo();
    expect(grok.shell.canRedo, false);
  }, { skipReason: nodeSkip, node: false });
});
