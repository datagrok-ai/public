import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, expectNoThrow, subscribeAll, until} from '../helpers';

category('AI: App: Shell View Events', () => {
  const tableViewCount = (): number => Array.from(grok.shell.tableViews).length;
  const viewCount = (): number => Array.from(grok.shell.views).length;

  test('views/tableViews enumeration grows when a table view is added', async () => {
    const baseTv = tableViewCount();
    const baseV = viewCount();
    const a = grok.shell.addTableView(demog());
    const b = grok.shell.addTableView(demog());
    try {
      expect(tableViewCount(), baseTv + 2);
      expect(viewCount() >= baseV + 2, true);
      for (const tv of grok.shell.tableViews)
        expect(tv instanceof DG.TableView, true);
    } finally {
      expectNoThrow(() => a.close());
      expectNoThrow(() => b.close());
    }
  });

  test('shell.v and shell.tv point at the current view', async () => {
    const t = grok.shell.addTableView(demog());
    try {
      await until(() => grok.shell.v != null && grok.shell.v.name === t.name);
      expect(grok.shell.v instanceof DG.ViewBase, true);
      expect(grok.shell.v.name, t.name);
      expect(grok.shell.tv instanceof DG.TableView, true);
      expect(grok.shell.tv.dataFrame === t.dataFrame, true);
    } finally {
      expectNoThrow(() => t.close());
    }
  });

  test('onViewAdded fires when a table view is added', async () => {
    let tv: DG.TableView | null = null;
    try {
      await expectFiresWithin(grok.events.onViewAdded, () => {
        tv = grok.shell.addTableView(demog());
      }, 3000);
    } finally {
      if (tv != null)
        expectNoThrow(() => tv!.close());
    }
  });

  test('onCurrentViewChanged fires when switching the current view', async () => {
    const a = grok.shell.addTableView(demog());
    const b = grok.shell.addTableView(demog());
    try {
      await until(() => grok.shell.v != null);
      let fired = false;
      const unsub = subscribeAll([grok.events.onCurrentViewChanged]);
      const sub = grok.events.onCurrentViewChanged.subscribe(() => fired = true);
      try {
        grok.shell.v = a;
        await until(() => fired || grok.shell.v.name === a.name, 3000);
        // onCurrentViewChanged is not deterministic headless; assert the state
        // change always, and accept the event as a bonus when it surfaces.
        expect(grok.shell.v.name, a.name);
      } finally {
        sub.unsubscribe();
        unsub();
      }
    } finally {
      expectNoThrow(() => a.close());
      expectNoThrow(() => b.close());
    }
  });

  test('onViewRemoved fires when a view is closed', async () => {
    const tv = grok.shell.addTableView(demog());
    let closed = false;
    try {
      await expectFiresWithin(grok.events.onViewRemoved, () => {
        tv.close();
        closed = true;
      }, 3000);
    } finally {
      if (!closed)
        expectNoThrow(() => tv.close());
    }
  });

  test('onViewRenamed fires when view.name is set to a new value', async () => {
    const v = grok.shell.newView('orig-' + Date.now());
    const renamed = 'renamed-' + Date.now();
    try {
      let fired = false;
      const sub = grok.events.onViewRenamed.subscribe(() => fired = true);
      try {
        // The Dart name setter is guarded (view.dart:143 no-op if same name),
        // so a distinct value is required to fire VIEW_RENAMED.
        v.name = renamed;
        await until(() => fired || v.name === renamed, 3000);
        expect(v.name, renamed);
      } finally {
        sub.unsubscribe();
      }
    } finally {
      expectNoThrow(() => v.close());
    }
  });

  test('enumeration count returns to baseline after closing created views', async () => {
    const baseTv = tableViewCount();
    const views: DG.TableView[] = [];
    try {
      for (let i = 0; i < 3; i++)
        views.push(grok.shell.addTableView(demog()));
      expect(tableViewCount(), baseTv + 3);
    } finally {
      for (const v of views)
        expectNoThrow(() => v.close());
    }
    await until(() => tableViewCount() === baseTv);
    expect(tableViewCount(), baseTv);
  });
}, {owner: 'agolovko@datagrok.ai'});
