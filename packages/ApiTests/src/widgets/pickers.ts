import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {waitUntil, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

async function cancelPick(pick: () => Promise<unknown>, title: string): Promise<void> {
  const picked = pick();
  try {
    await waitUntil(() => DG.Dialog.getOpenDialogs().some((d) => d.title === title),
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

async function openDialog(title: string): Promise<DG.Dialog> {
  await waitUntil(() => DG.Dialog.getOpenDialogs().some((d) => d.title === title),
    `the '${title}' dialog did not open`, 15000);
  return DG.Dialog.getOpenDialogs().find((d) => d.title === title)!;
}

/** The label of the Demo share's node, once the dialog's tree shows it. */
async function demoShare(dialog: DG.Dialog, selector: string): Promise<HTMLElement> {
  await waitUntil(() => dialog.root.querySelector(selector) !== null, 'the Demo share is not in the tree', 15000);
  return dialog.root.querySelector(selector) as HTMLElement;
}

function closeDialogs(title: string): void {
  for (const d of DG.Dialog.getOpenDialogs()) {
    if (d.title === title)
      d.close();
  }
}

function treeLabel(dialog: DG.Dialog, text: string): HTMLElement | undefined {
  return Array.from(dialog.root.querySelectorAll<HTMLElement>('.d4-tree-view-item-label, .d4-tree-view-group-label'))
    .find((e) => e.textContent === text);
}

/** Expands a tree group by its label, unless it is: the tree keeps the folders a previous dialog opened. */
function expand(label: Element): void {
  const expander = label.parentElement!.querySelector('.d4-tree-view-tri') as HTMLElement;
  if (!expander.classList.contains('d4-tree-view-tri-expanded'))
    expander.click();
}

category('Widgets: Table pickers', () => {
  test('pickFileFromFiles cancel resolves null', async () => {
    await cancelPick(() => ui.pickFileFromFiles(), 'Select a file');
  });

  test('pickFileFromFiles resolves the file selected on OK', async () => {
    const picked = ui.pickFileFromFiles();
    try {
      const dialog = await openDialog('Select a file');
      expand(await demoShare(dialog, '.d4-tree-view-group-label[data-namespace="System:"][data-name="DemoFiles"]'));
      await waitUntil(() => treeLabel(dialog, 'demog.csv') !== undefined,
        'demog.csv is not listed in the Demo share', 15000);
      treeLabel(dialog, 'demog.csv')!.click();
      (dialog.root.querySelector('[name="button-OK"]') as HTMLElement).click();
      const file = await picked;
      expect(file instanceof DG.FileInfo, true);
      expect(file!.fullPath, 'System:DemoFiles/demog.csv');
      expect(file!.isFile, true);
    } finally {
      closeDialogs('Select a file');
    }
  });

  test('pickFileFromFiles with multiple resolves every file selected on OK', async () => {
    const picked = ui.pickFileFromFiles({multiple: true});
    try {
      const dialog = await openDialog('Select a file');
      expand(await demoShare(dialog, '.d4-tree-view-group-label[data-namespace="System:"][data-name="DemoFiles"]'));
      await waitUntil(() => treeLabel(dialog, 'demog.csv') !== undefined && treeLabel(dialog, 'cars.csv') !== undefined,
        'demog.csv and cars.csv are not listed in the Demo share', 15000);
      treeLabel(dialog, 'demog.csv')!.click();
      treeLabel(dialog, 'cars.csv')!.dispatchEvent(new MouseEvent('click', {bubbles: true, ctrlKey: true}));
      (dialog.root.querySelector('[name="button-OK"]') as HTMLElement).click();
      const files = await picked;
      expect(Array.isArray(files), true);
      expectArray(files!.map((f) => f.fullPath).sort(), ['System:DemoFiles/cars.csv', 'System:DemoFiles/demog.csv']);
      expect(files!.every((f) => f instanceof DG.FileInfo && f.isFile), true);
    } finally {
      closeDialogs('Select a file');
    }
  });

  test('pickFileFromFiles with a connection browses that share alone', async () => {
    const demo = await grok.dapi.connections.filter('shortName = "DemoFiles"').first();
    const picked = ui.pickFileFromFiles({connection: demo});
    try {
      const dialog = await openDialog('Select a file');
      const share = await demoShare(dialog, '.d4-tree-view-group-label[data-name="DemoFiles"]');
      expect(dialog.root.querySelectorAll('.d4-tree-view-group-label[data-name]').length, 1,
        'only the connection\'s share');
      expand(share);
      await waitUntil(() => treeLabel(dialog, 'demog.csv') !== undefined,
        'demog.csv is not listed in the Demo share', 15000);
      treeLabel(dialog, 'demog.csv')!.click();
      (dialog.root.querySelector('[name="button-OK"]') as HTMLElement).click();
      const file = await picked;
      expect(file instanceof DG.FileInfo, true);
      expect(file!.fullPath, 'System:DemoFiles/demog.csv');
      expect(file!.path, 'demog.csv');
    } finally {
      closeDialogs('Select a file');
    }
  });

  test('pickFileFromFiles with a connection in folder mode lists folders and resolves the share\'s root', async () => {
    const demo = await grok.dapi.connections.filter('shortName = "DemoFiles"').first();
    const picked = ui.pickFileFromFiles({connection: demo, mode: 'folder'});
    try {
      const dialog = await openDialog('Select a directory');
      const share = await demoShare(dialog, '.d4-tree-view-group-label[data-name="DemoFiles"]');
      expand(share);
      await waitUntil(() => treeLabel(dialog, 'geo') !== undefined,
        'the geo folder is not listed in the Demo share', 15000);
      expect(treeLabel(dialog, 'demog.csv') === undefined, true, 'a file listed in folder mode');
      share.click();
      (dialog.root.querySelector('[name="button-OK"]') as HTMLElement).click();
      const folder = await picked;
      expect(folder!.isDirectory, true);
      expect(folder!.path, '');
      expect(folder!.fullPath, 'System:DemoFiles');
    } finally {
      closeDialogs('Select a directory');
    }
  });

  test('pickTableFromFiles cancel resolves null', async () => {
    await cancelPick(() => ui.pickTableFromFiles(), 'Select a file');
  });

  test('pickTableFromQuery cancel resolves null', async () => {
    await cancelPick(() => ui.pickTableFromQuery(), 'Select a database query');
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
      await waitUntil(() => cg.dfColumns.filter.trueCount === expected,
        'search did not filter the column list');
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
