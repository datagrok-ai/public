import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

// npm datagrok-api@1.27.9 typings predate the 1.28 ui.pickTableFrom* wrappers
const pickers = ui as any;

async function cancelPick(pick: () => Promise<DG.DataFrame | null>, title: string): Promise<void> {
  const picked = pick();
  try {
    await awaitCheck(() => DG.Dialog.getOpenDialogs().some((d) => d.title === title),
      `the '${title}' dialog did not open`, 15000);
    const dialog = DG.Dialog.getOpenDialogs().find((d) => d.title === title)!;
    (dialog.root.querySelector('[name="button-CANCEL"]') as HTMLElement).click();
    expect(await picked, null);
  } finally {
    for (const d of DG.Dialog.getOpenDialogs()) {
      if (d.title === title)
        d.close();
    }
  }
}

category('Widgets: Table pickers', () => {
  test('pickTableFromFiles cancel resolves null', async () => {
    await cancelPick(() => pickers.pickTableFromFiles(), 'Select a file');
  });

  test('pickTableFromQuery cancel resolves null', async () => {
    await cancelPick(() => pickers.pickTableFromQuery(), 'Select a database query');
  });
}, {owner: 'askalkin@datagrok.ai'});

category('Widgets: ColumnGrid', () => {
  let df: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(20);
  });

  test('popup constructs over all columns', async () => {
    const cg = DG.ColumnGrid.popup(df);
    try {
      expect(cg.root.classList.contains('d4-column-grid'), true);
      expect(cg.dfColumns.rowCount, df.columns.length);
    } finally {
      cg.close();
    }
  });

  test('current row change fires and currentColumn answers', async () => {
    const cg = DG.ColumnGrid.popup(df);
    try {
      let fired = 0;
      const sub = cg.dfColumns.onCurrentRowChanged.subscribe(() => fired++);
      try {
        cg.dfColumns.currentRowIdx = 2;
      } finally {
        sub.unsubscribe();
      }
      expect(fired, 1);
      expect(cg.currentColumn.name, cg.getCol(2).name);
    } finally {
      cg.close();
    }
  });

  test('popup filter callback receives a wrapped column', async () => {
    const seen: string[] = [];
    const cg = DG.ColumnGrid.popup(df, {filter: (c) => {
      seen.push(typeof c.name === 'string' ? c.name : '<unwrapped>');
      return c.isNumerical;
    }});
    try {
      expect(seen.length >= df.columns.length, true);
      expect(seen.every((name) => df.columns.contains(name)), true);
      const names: string[] = [];
      for (let i = 0; i < cg.dfColumns.rowCount; i++)
        names.push(cg.getCol(i).name);
      expectArray(names, df.columns.toList().filter((c) => c.isNumerical).map((c) => c.name));
    } finally {
      cg.close();
    }
  });

  test('filter setter re-filters with wrapped columns', async () => {
    const cg = DG.ColumnGrid.popup(df);
    try {
      cg.filter = (c) => c.name === 'age';
      cg.filterColumns();
      expect(cg.dfColumns.filter.trueCount, 1);
    } finally {
      cg.close();
    }
  });

  test('search box filters the list', async () => {
    const cg = DG.ColumnGrid.popup(df);
    const host = ui.div([cg.root]);
    document.body.append(host);
    try {
      cg.showSearch = true;
      const search = cg.root.querySelector('input') as HTMLInputElement;
      expect(search != null, true);
      search.value = 'age';
      search.dispatchEvent(new Event('input', {bubbles: true}));
      const expected = df.columns.names().filter((n) => n.toLowerCase().includes('age')).length;
      await awaitCheck(() => cg.dfColumns.filter.trueCount === expected,
        'search did not filter the column list', 3000);
    } finally {
      cg.close();
      host.remove();
    }
  });

  test('close detaches the root', async () => {
    const cg = DG.ColumnGrid.popup(df);
    const host = ui.div([cg.root]);
    document.body.append(host);
    try {
      expect(cg.root.isConnected, true);
      cg.close();
      expect(cg.root.isConnected, false);
    } finally {
      host.remove();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
