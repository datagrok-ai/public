import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray, awaitCheck} from '@datagrok-libraries/test/src/test';
import {openRunComparison} from '../components/RunComparison/comparison-views';
import {entryFromDataFrame} from '../components/RunComparison/entry-extraction';

const table = (name: string, values: number[]) => {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1, 2]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', values),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'same', [1, 2, 3]),
  ]);
  df.name = name;
  return df;
};

// the index pickers sit above the tabs; the Compare list, toggles and chart in the shown tab
const shownPane = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLElement>('[role="tabpanel"]'))
  .find((pane) => pane.style.display !== 'none') ?? null;
const indexSelects = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLElement>('.c2-comparison-table'))
  .filter((table) => !table.closest('[role="tabpanel"]'))
  .flatMap((table) => Array.from(table.querySelectorAll<HTMLElement>('.c2-comparison-row')))
  .map((row) => row.querySelector('select')!);
const targetRows = (root: HTMLElement) =>
  Array.from(shownPane(root)?.querySelectorAll<HTMLElement>('.c2-comparison-table .c2-comparison-row') ?? []);
const rowName = (row: HTMLElement) => row.querySelector('span[title]')!.getAttribute('title');
const targetNames = (root: HTMLElement) => targetRows(root).map(rowName);
const targetRow = (root: HTMLElement, name: string) => targetRows(root).find((row) => rowName(row) === name)!;
const selectedNames = (root: HTMLElement) => targetRows(root)
  .filter((row) => row.classList.contains('c2-comparison-row-selected')).map(rowName);
const chartElement = (root: HTMLElement) =>
  shownPane(root)?.querySelector('dg-viewer') as (HTMLElement & {dataFrame?: DG.DataFrame}) | null;
const chartFrame = (root: HTMLElement) => chartElement(root)?.dataFrame;
const hideEqual = (root: HTMLElement, value: boolean) => {
  const toggle = Array.from(shownPane(root)!.querySelectorAll<HTMLElement & {value: boolean}>('dg-toggle-input'))
    .find((el) => el.textContent?.includes('Hide equal values'))!;
  toggle.value = value;
};
const tabHeaders = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLElement>('[role="tab"]'));
const tabTitles = (root: HTMLElement) => tabHeaders(root).map((tab) => tab.getAttribute('name'));
const shownTitle = (root: HTMLElement) => shownPane(root)?.getAttribute('name');
const tab = (root: HTMLElement, title: string) =>
  tabHeaders(root).find((header) => header.getAttribute('name') === title)!;
// a hidden chart's canvas collapses; shown again, it must fill the chart element
const canvasFills = (chart: HTMLElement) => {
  const widths = Array.from(chart.querySelectorAll('canvas')).map((c) => c.getBoundingClientRect().width);
  const width = Math.max(0, ...widths);
  return width > 0 && Math.abs(width - chart.getBoundingClientRect().width) < 20;
};
// the lowest point of the element any clipping ancestor still shows
const visibleBottom = (el: HTMLElement) => {
  let bottom = el.getBoundingClientRect().bottom;
  for (let parent = el.parentElement; parent; parent = parent.parentElement) {
    if (getComputedStyle(parent).overflowY !== 'visible')
      bottom = Math.min(bottom, parent.getBoundingClientRect().bottom);
  }
  return bottom;
};
const newTab = (root: HTMLElement) => root.querySelector<HTMLElement>('[aria-label="New tab"]')!.click();

const pick = (select: HTMLSelectElement, value: string) => {
  select.value = value;
  select.dispatchEvent(new Event('change'));
};

// closes the view when setup fails: the caller's finally only covers a returned handle
async function openWithIndexes(...dfs: DG.DataFrame[]) {
  const handle = openRunComparison({entries: dfs.map(entryFromDataFrame)});
  const root = handle.view.root;
  try {
    await awaitCheck(() => indexSelects(root).length === dfs.length, 'index pickers not shown', 3000);
    indexSelects(root).forEach((select) => pick(select as HTMLSelectElement, 'time'));
    await awaitCheck(() => targetRows(root).length > 0, 'no values to compare after picking indexes', 3000);
  } catch (e) {
    handle.view.close();
    throw e;
  }
  return handle;
}

category('RunComparison: component', () => {
  test('values are listed once indexes are picked, and a click charts one', async () => {
    const handle = openRunComparison({entries: [table('A', [1, 2, 3]), table('B', [2, 3, 4])]
      .map(entryFromDataFrame)});
    const root = handle.view.root;
    try {
      await awaitCheck(() => indexSelects(root).length === 2, 'index pickers not shown', 3000);
      expect(targetRows(root).length, 0);
      indexSelects(root).forEach((select) => pick(select as HTMLSelectElement, 'time'));
      await awaitCheck(() => targetNames(root).includes('value'), 'value not listed', 3000);

      targetRow(root, 'value').click();
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: value', 'value not charted', 5000);
      expectArray(selectedNames(root), ['value']);
      expect(chartFrame(root)!.rowCount, 6);

      targetRow(root, 'value').click();
      await awaitCheck(() => chartElement(root) == null, 'chart not cleared', 3000);
    } finally {
      handle.view.close();
    }
  });

  test('equal values are hidden until the toggle is turned off', async () => {
    const handle = await openWithIndexes(table('A', [1, 2, 3]), table('B', [2, 3, 4]));
    const root = handle.view.root;
    try {
      expectArray(targetNames(root), ['value']);
      hideEqual(root, false);
      await awaitCheck(() => targetNames(root).length === 2, 'equal value not shown', 3000);
      expectArray(targetNames(root), ['same', 'value']);
    } finally {
      handle.view.close();
    }
  });

  test('shift click charts several values, anchor first', async () => {
    const handle = await openWithIndexes(table('A', [1, 2, 3]), table('B', [2, 3, 4]));
    const root = handle.view.root;
    try {
      hideEqual(root, false);
      await awaitCheck(() => targetNames(root).length === 2, 'equal value not shown', 3000);
      targetRow(root, 'value').click();
      targetRow(root, 'same').dispatchEvent(new MouseEvent('click', {bubbles: true, shiftKey: true}));
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: multiple values',
        'values not charted together', 5000);
      expectArray(selectedNames(root), ['same', 'value']);
      const valueColumns = chartFrame(root)!.columns.names().filter((name) => name === 'value' || name === 'same');
      expectArray(valueColumns, ['value', 'same']);
    } finally {
      handle.view.close();
    }
  });

  test('a removed run loses its picks', async () => {
    const first = table('A', [1, 2, 3]);
    const handle = await openWithIndexes(first, table('B', [2, 3, 4]), table('C', [3, 4, 5]));
    const root = handle.view.root;
    try {
      root.querySelector<HTMLElement>('[aria-label="Remove from comparison"]')!.click();
      await awaitCheck(() => indexSelects(root).length === 2, 'run not removed', 3000);
      await handle.service.add({entries: [entryFromDataFrame(first)]});
      await awaitCheck(() => indexSelects(root).length === 3, 'run not added back', 3000);
      expectArray(indexSelects(root).map((select) => (select as HTMLSelectElement).value), ['time', 'time', '']);
    } finally {
      handle.view.close();
    }
  });

  test('export opens a copy of the chart data in the workspace', async () => {
    const handle = await openWithIndexes(table('A', [1, 2, 3]), table('B', [2, 3, 4]));
    const root = handle.view.root;
    let exported: DG.TableView | undefined;
    try {
      targetRow(root, 'value').click();
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: value', 'value not charted', 5000);
      Array.from(shownPane(root)!.querySelectorAll<HTMLElement>('.c2-comparison-link'))
        .find((link) => link.textContent?.includes('Export'))!.click();
      const openButton = () => Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog button'))
        .find((button) => button.textContent?.trim() === 'OPEN IN WORKSPACE');
      await awaitCheck(() => openButton() != null, 'export dialog not shown', 3000);
      openButton()!.click();
      await awaitCheck(() => grok.shell.tv?.dataFrame?.name === 'Comparison: value', 'export not opened', 5000);
      exported = grok.shell.tv;
      expect(exported.dataFrame !== chartFrame(root), true);
      expect(exported.dataFrame.rowCount, 6);
    } finally {
      document.querySelectorAll('.d4-dialog').forEach((el) => el.remove());
      if (exported) {
        exported.close();
        grok.shell.closeTable(exported.dataFrame);
      }
      handle.view.close();
    }
  });

  test('tabs keep their own selection and toggles over shared runs and picks', async () => {
    const handle = await openWithIndexes(table('A', [1, 2, 3]), table('B', [2, 3, 4]));
    const root = handle.view.root;
    try {
      expectArray(tabTitles(root), ['Comparison 1']);
      targetRow(root, 'value').click();
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: value', 'value not charted', 5000);
      await awaitCheck(() => tabTitles(root)[0] === 'value', 'tab not named after the charted value', 3000);
      const firstChart = chartElement(root);

      newTab(root);
      await awaitCheck(() => shownTitle(root) === 'Comparison 2', 'new tab not shown', 3000);
      expectArray(tabTitles(root), ['value', 'Comparison 2']);
      expectArray(selectedNames(root), []);
      expect(chartElement(root) == null, true);
      hideEqual(root, false);
      await awaitCheck(() => targetNames(root).length === 2, 'equal value not shown in the new tab', 3000);
      targetRow(root, 'same').click();
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: same', 'same not charted', 5000);
      await awaitCheck(() => tabTitles(root)[1] === 'same', 'second tab not named', 3000);

      tab(root, 'value').click();
      await awaitCheck(() => shownTitle(root) === 'value', 'first tab not shown', 3000);
      expectArray(targetNames(root), ['value']);
      expectArray(selectedNames(root), ['value']);
      expect(chartElement(root) === firstChart, true);
      expect(chartFrame(root)!.name, 'Comparison: value');
      await awaitCheck(() => canvasFills(firstChart!), 'the chart did not redraw at full size', 3000);

      // an index pick is shared: unsetting it leaves both tabs without column values
      indexSelects(root).forEach((select) => pick(select as HTMLSelectElement, ''));
      await awaitCheck(() => targetRows(root).length === 0, 'first tab still lists values', 3000);
      tab(root, 'Comparison 2').click();
      await awaitCheck(() => shownTitle(root) === 'Comparison 2', 'second tab not shown', 3000);
      expectArray(targetNames(root), []);
    } finally {
      handle.view.close();
    }
  });

  test('closing the shown tab shows its neighbor, and the last tab cannot be closed', async () => {
    const handle = await openWithIndexes(table('A', [1, 2, 3]), table('B', [2, 3, 4]));
    const root = handle.view.root;
    const close = (title: string) => tab(root, title).querySelector<HTMLElement>('[aria-label="Close"]')!.click();
    try {
      newTab(root);
      newTab(root);
      await awaitCheck(() => tabTitles(root).length === 3, 'tabs not added', 3000);
      tab(root, 'Comparison 2').click();
      await awaitCheck(() => shownTitle(root) === 'Comparison 2', 'second tab not shown', 3000);

      close('Comparison 2');
      await awaitCheck(() => tabTitles(root).length === 2, 'tab not closed', 3000);
      expect(shownTitle(root), 'Comparison 3');
      close('Comparison 3');
      await awaitCheck(() => tabTitles(root).length === 1, 'tab not closed', 3000);
      expect(shownTitle(root), 'Comparison 1');
      expect(tab(root, 'Comparison 1').querySelector('[aria-label="Close"]') == null, true);

      newTab(root);
      await awaitCheck(() => shownTitle(root) === 'Comparison 4', 'numbers of closed tabs are reused', 3000);
    } finally {
      handle.view.close();
    }
  });

  test('a tab taller than the view scrolls to its chart bottom', async () => {
    const values = Array.from({length: 8}, (_, i) => `v${i + 1}`);
    const tall = (name: string, shift: number) => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1, 2]),
        ...values.map((value, i) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, value, [i, i + 1, i + shift])),
      ]);
      df.name = name;
      return df;
    };
    const handle = await openWithIndexes(tall('A', 2), tall('B', 3));
    const root = handle.view.root;
    try {
      targetRow(root, 'v1').click();
      for (const value of values.slice(1)) {
        targetRow(root, value).dispatchEvent(new MouseEvent('click', {bubbles: true, shiftKey: true}));
        await awaitCheck(() => selectedNames(root).includes(value), `${value} not added`, 3000);
      }
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: multiple values',
        'values not charted together', 5000);
      const chart = chartElement(root)!;
      expect(chart.getBoundingClientRect().height > root.getBoundingClientRect().height, true);
      // only overflow auto/scroll scrolls with the wheel; overflow hidden can scroll from code only
      for (let el = chart.parentElement; el; el = el.parentElement) {
        if (['auto', 'scroll'].includes(getComputedStyle(el).overflowY))
          el.scrollTop = el.scrollHeight;
      }
      const {bottom} = chart.getBoundingClientRect();
      expect(Math.round(visibleBottom(chart)) >= Math.round(bottom) - 1, true);
    } finally {
      handle.view.close();
    }
  });

  test('tabs keep their state while fewer than two runs are compared', async () => {
    // one table name, like runs of one model: a column target's key includes the table name
    const handle = await openWithIndexes(table('Data', [1, 2, 3]), table('Data', [2, 3, 4]));
    const root = handle.view.root;
    try {
      targetRow(root, 'value').click();
      await awaitCheck(() => tabTitles(root)[0] === 'value', 'tab not named after the charted value', 3000);
      newTab(root);
      await awaitCheck(() => shownTitle(root) === 'Comparison 2', 'new tab not shown', 3000);
      hideEqual(root, false);
      await awaitCheck(() => targetNames(root).length === 2, 'equal value not shown', 3000);

      root.querySelector<HTMLElement>('[aria-label="Remove from comparison"]')!.click();
      await awaitCheck(() => root.textContent!.includes('Add at least two runs'), 'one run left, no message', 3000);
      await handle.service.add({entries: [entryFromDataFrame(table('Data', [5, 6, 7]))]});
      await awaitCheck(() => indexSelects(root).length === 2, 'runs not shown again', 3000);
      indexSelects(root).forEach((select) => pick(select as HTMLSelectElement, 'time'));

      await awaitCheck(() => targetNames(root).length === 2, 'the second tab lost Hide equal values off', 3000);
      await awaitCheck(() => tabTitles(root)[0] === 'value', 'the first tab lost its value', 3000);
      tab(root, 'value').click();
      await awaitCheck(() => shownTitle(root) === 'value', 'the first tab lost its value', 3000);
      expectArray(selectedNames(root), ['value']);
      await awaitCheck(() => chartFrame(root)?.name === 'Comparison: value', 'value not charted again', 5000);
    } finally {
      handle.view.close();
    }
  });
});

