import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';

// Source scenario: General/table-manager.md (automatable slice).
//
// The original scenario's canvas-grid interactions — Shift+click multi-select of
// table rows and the resulting "N tables" context submenu, plus the "Show > All"
// attribute-column toggle and "Open as table" right-click action — are driven by
// pixel coordinates on the Table Manager's d4 grid and have no stable DOM/API
// hook, so they live in `table-manager-ui.md`.
//
// What IS robustly automatable: opening the Table Manager (View | Tables /
// Alt+T -> `grok.shell.windows.showTables`) and verifying it lists every open
// table. Verified live on dev 2026-06-16: the pane renders as `[name="Tables"]`
// containing a `.grok-tables-manager` grid whose rows expose each table name in
// the DOM.

test.use(specTestOptions);

const TABLE_NAMES = ['ttManagerA', 'ttManagerB', 'ttManagerC'];

test('Table Manager lists open tables', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  // Precondition: load several datasets (scenario asks for 3-4).
  await page.evaluate(() => (window as any).grok.shell.closeAll());
  await page.waitForFunction(() => (window as any).grok.shell.tables.length === 0, null, {timeout: 15_000});
  const opened = await page.evaluate((names) => {
    const grok = (window as any).grok;
    for (const n of names) {
      const t = grok.data.testData('demog', 40);
      t.name = n;
      grok.shell.addTableView(t);
    }
    return {count: grok.shell.tables.length, names: grok.shell.tableNames};
  }, TABLE_NAMES);

  expect(opened.count).toBe(TABLE_NAMES.length);
  for (const n of TABLE_NAMES) expect(opened.names).toContain(n);

  // Open the Table Manager (View | Tables / Alt+T).
  await page.evaluate(() => (window as any).grok.shell.windows.showTables = true);
  await page.locator('[name="Tables"]').waitFor({state: 'attached', timeout: 15_000});
  const manager = await page.evaluate((names) => {
    const grok = (window as any).grok;
    const textVisible = (txt: string) =>
      [...document.querySelectorAll('div,span,td')]
        .some((e) => e.textContent === txt && (e as HTMLElement).offsetParent !== null);
    return {
      getter: grok.shell.windows.showTables,
      panePresent: !!document.querySelector('[name="Tables"]'),
      gridPresent: !!document.querySelector('.grok-tables-manager'),
      namesRendered: names.map((n: string) => ({n, found: textVisible(n)})),
    };
  }, TABLE_NAMES);

  expect(manager.getter).toBe(true);
  expect(manager.panePresent).toBe(true);
  expect(manager.gridPresent).toBe(true);
  for (const row of manager.namesRendered)
    expect(row.found, `table "${row.n}" should be listed in the Table Manager`).toBe(true);

  // Toggling View | Tables off hides the pane (the command is a simple toggle).
  await softStep('View | Tables toggles the manager off', async () => {
    await page.evaluate(() => (window as any).grok.shell.windows.showTables = false);
    await page.waitForFunction(() => (window as any).grok.shell.windows.showTables === false,
      null, {timeout: 15_000});
  });

  // Cleanup.
  await page.evaluate(() => (window as any).grok.shell.closeAll());
});
