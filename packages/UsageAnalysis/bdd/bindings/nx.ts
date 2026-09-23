/* The steps the NX chain needs (features/viewers/nx): the Link Tables dialog's tables and key
   pairs, the line chart's split selectors, the projects the chain saves through the ribbon's Save
   dialog and reopens by their friendly name, the Formula Lines dialog, the Scaffold Tree filter
   card, layouts kept by name, and the claims that compare what a viewer shows, or how many rows
   pass, with the table itself, across views and across a project round trip. */
import {Locator, Page} from '@playwright/test';
import {dataset, element, Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, type ElementRef, expect, gestures, pollMs, projectSaveWindow, viewers} from '@datagrok-libraries/bdd/runtime';

dataset('spgi-3624', {path: 'System:DemoFiles/chem/SPGI.csv', description: 'the full SPGI demo table: 3624 molecules, 88 columns, opened as "SPGI"'});

// the Add New Column dialog's name field has neither a label nor a placeholder
element('new column name input', {selector: '[name="input-Add-New-Column---Name"]'});
// the column list's own checkbox beside its search box: it checks or unchecks every column at once
element('all columns checkbox', {selector: '[name="dialog-Order-or-Hide-Columns"] .d4-column-grid > .ui-input-root input[type="checkbox"]'});

const linkDialog = (page: Page): Locator => page.locator('[name="dialog-Link-Tables"]').filter({visible: true}).last();

/** A pair of key columns of the Link Tables dialog: row N holds the column of the left table and
 * the column of the right one, each a Dart column selector. */
export const setKeyPair = When('user sets key columns {int} of the Link Tables dialog to {string} and {string}',
  async (page: Page, row: number, left: string, right: string) => {
    for (const [side, column] of [[0, left], [1, right]] as [number, string][]) {
      const selector = linkDialog(page).locator(`[name="div-selectKeyCol${side}Row${row}"] .d4-column-selector`).first();
      await expect(selector, `key column ${side === 0 ? 'left' : 'right'} of row ${row} in the Link Tables dialog`).toBeVisible({timeout: 5000});
      await gestures.openColumnSelector(page, selector);
      await gestures.pickInColumnGrid(page, column, `key column row ${row}`, selector);
      await expect(selector.locator('.d4-column-selector-column')).toHaveText(column, {timeout: 5000});
    }
  }, {tier: 'ui', description: 'picks the left and the right column of the key row N (1-based; "Add" makes a row) through their column pickers and reads each back'});

/** The two table choices of a link: the left one is labelled "Tables", the right one has no label. */
export const setLinkTables = When('user sets the tables of the Link Tables dialog to {string} and {string}',
  async (page: Page, left: string, right: string) => {
    for (const [side, table] of [['Left', left], ['Right', right]]) {
      const select = linkDialog(page).locator(`select[name="input-selectTable${side}"]`).first();
      await select.selectOption({label: table});
      await expect(select, `the ${side.toLowerCase()} table of the Link Tables dialog`).toHaveValue(table, {timeout: 5000});
    }
  }, {tier: 'ui', description: 'chooses the left and the right table of the link being made and reads both back'});

/** The line chart's split selectors sit in its top-left corner, revealed on hover: the ones that
 * hold a split, and the empty "add-split" one after them that takes the next. */
export const addLineChartSplit = When('user adds {string} to the splits of {widget}', async (page: Page, column: string, target: ElementRef) => {
  const root = await viewers.viewerLocator(page, target);
  await root.hover();
  const add = root.locator('[name="add-split"]').first();
  await expect(add, 'the empty split selector of the line chart').toBeVisible({timeout: 5000});
  await viewers.snapshot(page, target);
  await gestures.openColumnSelector(page, add);
  await gestures.pickInColumnGrid(page, column, 'the split selector of the line chart');
  await viewers.settle(page, target);
}, {tier: 'ui', description: 'hovers the chart, opens its empty split selector and picks the column in the popup'});

/** SPGI.csv is 11 MB: the view it opens as arrives well past a check's default wait. */
export const tableViewOpened = Then('the {string} table view should open with {int} rows', async (page: Page, name: string, rows: number) => {
  await expect.poll(() => page.evaluate((n) => {
    const tv = Array.from((window as any).grok.shell.tableViews as any[]).find((v: any) => v.dataFrame?.name === n);
    return tv ? tv.dataFrame.rowCount : -1;
  }, name), {timeout: pollMs(120000), message: `the rows of the "${name}" table view (-1: not open)`}).toBe(rows);
  await viewers.settleAll(page);
}, {description: 'polls for a table view of that table for as long as a large file takes to load (up to two minutes), then its row count'});

/** Opened through the project API: the Dashboards gallery search does not find a name that holds
 * "-", and a run-suffixed name holds several. */
export const openSavedProject = When('user opens the project saved as {string}', async (page: Page, name: string) => {
  projectSaveWindow(page, false);
  await page.evaluate(async (n) => {
    const found = await (window as any).grok.dapi.projects.filter(`friendlyName = "${n}"`).list();
    if (found.length !== 1)
      throw new Error(`${found.length} projects are named "${n}" on the server`);
    const project = await (window as any).grok.dapi.projects.find(found[0].id);
    await project.open();
  }, name);
  await page.waitForFunction(() => (window as any).grok.shell.tv?.dataFrame != null, null, {timeout: pollMs(120000)});
  await viewers.settleAll(page);
}, {tier: 'api', description: 'ends the window a save opened for its preview\'s console noise; the project whose friendly name is this (exactly one must exist); a project with data sync reads its files again, so this waits up to two minutes for a table view, then for its viewers to settle'});

export const projectLinks = Then('the project saved as {string} should link {string}',
  async (page: Page, name: string, links: string) => {
    const want = links.split(/\s*;\s*/).filter(Boolean).sort();
    await expect.poll(() => page.evaluate(async (n) => {
      const found = await (window as any).grok.dapi.projects.filter(`friendlyName = "${n}"`).list();
      if (found.length !== 1)
        return [`${found.length} projects named "${n}"`];
      const project = await (window as any).grok.dapi.projects.find(found[0].id);
      return ((project.options?.['table links'] ?? []) as any[])
        .map((l) => `${l.table1Name} -> ${l.table2Name} by ${l.keyColumns1.join(', ')} = ${l.keyColumns2.join(', ')} as ${l.linkTypes.join(', ')}`).sort();
    }, name), {message: `the table links the project "${name}" holds on the server`}).toEqual(want);
  }, {description: 'the links stored in the saved project, "A -> B by keys = keys as type", ";"-separated in any order'});

/** The rows of a table that pass its filter (what a link brought) and a viewer formula filter
 * of the "<text column> is <value> and <numeric column> is below <n>" shape — what a viewer
 * bound to that table with that formula should draw. */
async function passingRows(page: Page, table: string, textCol: string, value: string, numCol: string, below: number): Promise<number> {
  return page.evaluate(([t, tc, v, nc, b]) => {
    const df = (window as any).grok.shell.tables.find((x: any) => x.name === t);
    if (!df)
      throw new Error(`no table "${t}" is open`);
    const text = df.col(tc);
    const num = df.col(nc);
    if (!text || !num)
      throw new Error(`table "${t}" has no "${text ? nc : tc}" column`);
    if (!text.categories.includes(v))
      throw new Error(`"${tc}" holds no "${v}"; it holds: ${text.categories.join(', ')}`);
    let n = 0;
    for (let i = 0; i < df.rowCount; i++)
      if (df.filter.get(i) && text.get(i) === v && !num.isNone(i) && num.get(i) < b)
        n++;
    return n;
  }, [table, textCol, value, numCol, below] as [string, string, string, string, number]);
}

export const viewerShowsFormulaRows = Then('{widget} should show the rows of table {string} that pass the filter where {string} is {string} and {string} is below {float}',
  async (page: Page, target: ElementRef, table: string, textCol: string, value: string, numCol: string, below: number) => {
    let want = -1;
    await expect.poll(async () => {
      want = await passingRows(page, table, textCol, value, numCol, below);
      return (await viewers.readValue(page, target, 'rows shown')) === want;
    }, {message: `the "rows shown" reading against the rows of "${table}"`}).toBe(true).catch(async () => {
      throw new Error(`the viewer shows ${await viewers.readValue(page, target, 'rows shown')} rows; ${want} rows of "${table}" pass its filter and the formula`);
    });
  }, {description: 'the viewer\'s "rows shown" equals the rows of the table that pass the table filter (what the links put there) and hold the value and a number below the bound'});

/** Every table row of the Save project dialog carries its own Data sync switch; a row scrolled out
 * of the dialog's list is still read. */
export const saveDialogDataSync = Then('the Save project dialog should save the tables {string} with data sync',
  async (page: Page, tables: string) => {
    const want = tables.split(/\s*,\s*/).filter(Boolean).map((t) => `${t}: on`).sort();
    await expect.poll(() => page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('[name="dialog-Save-project"]')).pop();
      return Array.from(dialog?.querySelectorAll('.grok-project-move-entity-row') ?? [])
        .filter((row) => row.querySelector('[name="icon-table"]'))
        .map((row) => `${row.querySelector('label')?.textContent?.trim()}: ${row.querySelector('[name="input-host-Data-sync"] [role="switch"]')?.getAttribute('aria-checked') === 'true' ? 'on' : 'off'}`);
    }).then((rows) => rows.sort()), {message: 'the tables of the Save project dialog and their Data sync switches'}).toEqual(want);
  }, {description: 'the dialog lists exactly these tables (comma-separated), each with its Data sync switch on'});

/** The tab of the current view: its handle is named after the view, and a right-click on it opens
 * the view's own menu (View, Table, Dashboard). */
export const pickViewTabMenu = When('user picks {string} from the context menu of the current view tab', async (page: Page, path: string) => {
  const name = await page.evaluate(() => String((window as any).grok.shell.v?.name ?? ''));
  const tab = page.locator(`[name="view-handle: ${name}"]`).filter({visible: true}).last();
  await expect(tab, `the tab of the current view "${name}"`).toBeVisible({timeout: 5000});
  await tab.click({button: 'right'});
  await viewers.pickMenuPath(page, path);
}, {tier: 'ui', description: 'right-clicks the tab of the view in front and picks the path in the menu it opens'});

/** A project with several views of one table: the md goes to the last of them. */
export const switchToLastView = When('user switches to the last table view of {string}', async (page: Page, table: string) => {
  const name = await page.evaluate((t) => {
    const views = Array.from((window as any).grok.shell.tableViews as any[]).filter((v) => v.dataFrame?.name === t);
    if (views.length === 0)
      throw new Error(`no table view of "${t}" is open`);
    const last = views[views.length - 1];
    (window as any).grok.shell.v = last;
    return String(last.name);
  }, table);
  await page.waitForFunction((n) => String((window as any).grok.shell.v?.name) === n, name);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'the view of that table opened last, made current'});

/** SPGI has ninety columns and the view shows a dozen: a header the md right-clicks is scrolled
 * into view first, as a user drags the scroll bar to it. */
export const scrollGridTo = When('user scrolls the grid to the {string} column', async (page: Page, column: string) => {
  await page.evaluate((c) => {
    const grid = (window as any).grok.shell.tv?.grid;
    if (!grid?.dataFrame.col(c))
      throw new Error(`the grid of the current view has no "${c}" column`);
    grid.scrollToCell(c, 0);
  }, column);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'grid.scrollToCell on the first row of the column, then the viewers settle'});

/** "Not freezing, cannot be broken": every viewer of the view in front has finished its render and
 * reports no error of its own. */
export const noViewerError = Then('no viewer of the current view should report an error', async (page: Page) => {
  await viewers.settleAll(page);
  let seen: string[] = [];
  await expect.poll(async () => (seen = await page.evaluate(() => {
    const root = (window as any).grok.shell.v?.root as HTMLElement | undefined;
    const out: string[] = [];
    for (const el of Array.from(root?.querySelectorAll('[name^="viewer-"]') ?? [])) {
      if ((el as HTMLElement).getBoundingClientRect().width === 0 || el.parentElement?.closest('[name^="viewer-"]'))
        continue;
      let status: any;
      try {
        status = (window as any).__bdd.viewerOf(el)?.getWidgetStatus?.();
      }
      catch {
        status = undefined;
      }
      out.push(`${el.getAttribute('name')}: ${status ? String(status.error ?? '') : 'reports no status'}`);
    }
    return out;
  })).filter((s) => !s.endsWith(': ')), {timeout: pollMs(5000), message: 'viewers of the current view reporting an error'}).toEqual([]);
  expect(seen.length, 'the viewers of the current view').toBeGreaterThan(0);
}, {description: 'settles every viewer, then reads the error of each top-level viewer drawn in the current view; a viewer that reports no status fails, and at least one viewer must be there'});

export const showsFilteredRows = Then('{widget} should show every row that passes the filter of its table', async (page: Page, target: ElementRef) => {
  let got: {shown: number; passing: number; table: string} = {shown: -1, passing: -1, table: ''};
  await expect.poll(async () => {
    got = await viewers.onViewer(page, target, (el) => {
      const w = (window as any).__bdd.viewerOf(el);
      const df = w.dataFrame;
      // a scatter plot counts the markers it drew: a row with no value on an axis, or none above
      // zero on a logarithmic one, has no place on the plot
      const axes = ['x', 'y'].filter((a) => w.type === 'Scatter plot' && df.col(w.props[`${a}ColumnName`]))
        .map((a) => ({col: df.col(w.props[`${a}ColumnName`]), log: w.props[`${a}AxisType`] === 'logarithmic'}));
      let passing = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i) && axes.every((a) => !a.col.isNone(i) && (!a.log || a.col.get(i) > 0)))
          passing++;
      return {shown: Number(w.getWidgetStatus()?.values?.['rows shown'] ?? -1), passing, table: String(df.name)};
    });
    return got.shown === got.passing;
  }, {message: `the "rows shown" reading of ${target.phrase} against the rows its table lets through`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} shows ${got.shown} rows; ${got.passing} rows of "${got.table}" pass its filter`);
  });
}, {description: 'the viewer\'s "rows shown" equals the filter count of the table it is bound to (what the links left)'});

/** Pack and zoom by filter: the axes hug the rows that pass the filter. */
export const zoomedToFilter = Then('{widget} should be zoomed to the rows that pass the filter', async (page: Page, target: ElementRef) => {
  let why = '';
  await expect.poll(async () => (why = await viewers.onViewer(page, target, (el) => {
    const w = (window as any).__bdd.viewerOf(el);
    const df = w.dataFrame;
    const values = w.getWidgetStatus()?.values ?? {};
    const problems: string[] = [];
    const wide: string[] = [];
    for (const axis of ['x', 'y']) {
      const col = df.col(w.props[`${axis}ColumnName`]);
      const log = w.props[`${axis}AxisType`] === 'logarithmic';
      let fmin = Infinity; let fmax = -Infinity; let amin = Infinity; let amax = -Infinity;
      for (let i = 0; i < df.rowCount; i++) {
        if (col.isNone(i) || (log && col.get(i) <= 0))
          continue;
        const x = col.get(i);
        amin = Math.min(amin, x); amax = Math.max(amax, x);
        if (df.filter.get(i)) {
          fmin = Math.min(fmin, x); fmax = Math.max(fmax, x);
        }
      }
      const lo = Number(values[`${axis} axis min`]); const hi = Number(values[`${axis} axis max`]);
      const eps = (amax - amin) * 1e-6;
      if (!(lo <= fmin + eps && hi >= fmax - eps))
        problems.push(`${axis} axis ${lo}..${hi} leaves out the filtered ${fmin}..${fmax}`);
      if (fmax - fmin >= (amax - amin) * 0.9)
        wide.push(`${axis}: the filtered ${fmin}..${fmax} spans most of ${amin}..${amax}`);
      else if (!(hi - lo < amax - amin))
        problems.push(`${axis} axis ${lo}..${hi} spans the whole column ${amin}..${amax}, not the filtered ${fmin}..${fmax}`);
    }
    if (wide.length === 2)
      problems.push(`the data do not tell a zoom from none (${wide.join('; ')})`);
    return problems.join('; ');
  })) === '', {message: `the axes of ${target.phrase} against its filtered rows`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} is not zoomed to its filtered rows: ${why}`);
  });
}, {description: 'on each axis the range covers every value of the filtered rows (those above zero on a logarithmic axis), and where those span less than 90% of the column it is narrower than the whole column; filtered rows spanning most of the column on both axes fail, since they cannot tell a zoom from none'});

export const showsSelectedRows = Then('{widget} should show the selected rows of its table', async (page: Page, target: ElementRef) => {
  let got: {shown: number; selected: number; table: string} = {shown: -1, selected: -1, table: ''};
  await expect.poll(async () => {
    got = await viewers.onViewer(page, target, (el) => {
      const w = (window as any).__bdd.viewerOf(el);
      return {shown: Number(w.getWidgetStatus()?.values?.['rows shown'] ?? -1), selected: w.dataFrame.selection.trueCount, table: String(w.dataFrame.name)};
    });
    return got.shown === got.selected && got.selected > 0;
  }, {message: `the "rows shown" reading of ${target.phrase} against the selected rows of its table`}).toBe(true).catch(() => {
    throw new Error(`${target.phrase} shows ${got.shown} rows; "${got.table}" has ${got.selected} selected`);
  });
}, {description: 'a viewer whose Row Source is Selected: "rows shown" equals the table\'s selection, which must not be empty'});

export const clickViewTab = When('user clicks on the tab of the {string} view', async (page: Page, name: string) => {
  const tab = page.locator(`[name="view-handle: ${name}"]`).filter({visible: true}).last();
  await expect(tab, `the tab of the "${name}" view`).toBeVisible({timeout: 5000});
  await tab.click();
  await page.waitForFunction((n) => String((window as any).grok.shell.v?.name) === n, name, {timeout: pollMs(5000)});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the view\'s tab in the tab strip; done when that view is in front and its viewers have settled'});

// --- the Formula Lines dialog (PowerPack) -----------------------------------------------------------

const formulaDialog = (page: Page): Locator => page.locator('[name="dialog-Formula-Lines"]').filter({visible: true}).last();

/** The formulas the dialog's list holds on its current tab, read from the list's own grid. */
async function listedFormulas(page: Page): Promise<string[]> {
  return formulaDialog(page).locator('[name="viewer-Grid"]').filter({visible: true}).first().evaluate((el) => {
    const values = (window as any).__bdd.viewerOf(el).getWidgetStatus()?.values ?? {};
    return Object.keys(values).filter((k) => /^text of cell \d+ of formula$/.test(k))
      .sort((a, b) => Number(/\d+/.exec(a)![0]) - Number(/\d+/.exec(b)![0])).map((k) => String(values[k]));
  });
}

/** "Add new" > Line puts a line with a default formula in the list and selects it; its formula is
 * then typed into the editor's Line pane, which the list takes on blur. */
export const addFormulaLine = When('user adds the formula line {string} in the Formula Lines dialog', async (page: Page, formula: string) => {
  await viewers.settleAll(page);
  const dialog = formulaDialog(page);
  await dialog.locator('[name="button-Add-new"]').click();
  await viewers.pickMenuPath(page, 'Line');
  const editor = dialog.locator('[name="pane-Line"] textarea').first();
  await expect(editor, 'the formula editor of the new line').toBeVisible({timeout: 5000});
  await editor.fill(formula);
  await editor.press('Tab');
  await expect.poll(() => listedFormulas(page), {message: 'the formulas the Formula Lines dialog lists'}).toContain(formula);
}, {tier: 'ui', description: 'Add new > Line, the formula typed into the Line pane of the new line; done when the list shows it'});

export const setFormulaLineRange = When('user sets the range of the selected formula line to {string} .. {string}', async (page: Page, min: string, max: string) => {
  const host = formulaDialog(page).locator('[name="input-host-Range"]').first();
  const inputs = [host.locator('input').first(), host.locator('xpath=following-sibling::div[1]//input').first()];
  for (const [input, value] of [[inputs[0], min], [inputs[1], max]] as [Locator, string][]) {
    await input.fill(value);
    await input.press('Tab');
    await expect(input).toHaveValue(value, {timeout: 5000});
  }
}, {tier: 'ui', description: 'the min and the max field of the Range row of the Format pane'});

/** A line of the list is selected by a click on its formula; its formula is then retyped in the
 * editor's Line pane. */
export const editFormulaLine = When('user changes formula line {int} in the Formula Lines dialog to {string}', async (page: Page, row: number, formula: string) => {
  const dialog = formulaDialog(page);
  const grid = dialog.locator('[name="viewer-Grid"]').filter({visible: true}).first();
  const box = await grid.evaluate((el, r) => {
    const area = (window as any).__bdd.viewerOf(el).getWidgetStatus()?.hitAreas?.[`cell ${r} of formula`];
    const c = ((window as any).__bdd.viewerOf(el).getWidgetStatus()?.parts?.canvas ?? el.querySelector('canvas')).getBoundingClientRect();
    return area ? {x: c.x + area.x + area.width / 2, y: c.y + area.y + area.height / 2} : null;
  }, row);
  if (!box)
    throw new Error(`the Formula Lines dialog lists no line ${row}`);
  await page.mouse.click(box.x, box.y);
  const editor = dialog.locator('[name="pane-Line"] textarea').first();
  await expect(editor, 'the formula editor of the selected line').toBeVisible({timeout: 5000});
  await editor.fill(formula);
  await editor.press('Tab');
  await expect.poll(async () => (await listedFormulas(page))[row - 1], {message: `formula line ${row} of the Formula Lines dialog`}).toBe(formula);
}, {tier: 'ui', description: 'a click on the formula of line N in the list (the list puts the newest line first), the new formula typed into the Line pane; done when the list shows it'});

export const tableTagContains = Then('the {string} tag of the table should contain {string}', async (page: Page, tag: string, text: string) => {
  await expect.poll(() => page.evaluate((t) => String((window as any).grok.shell.tv?.dataFrame?.getTag(t) ?? ''), tag),
    {message: `the "${tag}" tag of the current table`}).toContain(text);
}, {description: 'the tag of the table in front — a dataframe formula line is kept in ".formula-lines"'});

// --- layouts the chain keeps by name ------------------------------------------------------------------

/** The chain applies a layout one scenario saved in a later one, on another view: the layouts stay
 * in the page, which the journey keeps. */
export const saveNamedLayout = When('user saves the layout of the current table view as {string}', async (page: Page, name: string) => {
  await viewers.settleAll(page);
  await page.evaluate((n) => {
    const w = window as any;
    w.__nxLayouts = {...(w.__nxLayouts ?? {}), [n]: w.grok.shell.tv.saveLayout()};
  }, name);
}, {tier: 'api', description: 'tv.saveLayout() kept under a name for "applies the layout" in this or a later scenario'});

export const applyNamedLayout = When('user applies the layout {string} to the current table view', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const layout = ((window as any).__nxLayouts ?? {})[n];
    if (!layout)
      throw new Error(`no layout "${n}" was saved; saved: ${Object.keys((window as any).__nxLayouts ?? {}).join(', ')}`);
    (window as any).grok.shell.tv.loadLayout(layout);
  }, name);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'tv.loadLayout of a layout saved by name, then every viewer settles'});

// --- the filter count, remembered across views ------------------------------------------------------

const passing = (page: Page): Promise<number> => page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.trueCount);

export const rememberPassing = When('user remembers how many rows pass the filter', async (page: Page) => {
  await viewers.settleAll(page);
  await page.evaluate(() => { (window as any).__nxPassing = (window as any).grok.shell.tv.dataFrame.filter.trueCount; });
}, {tier: 'api', description: 'the filter count of the table in front, for a claim in another view or after a change'});

export const passingAsRemembered = Then('as many rows as remembered should pass the filter', async (page: Page) => {
  const want = await page.evaluate(() => (window as any).__nxPassing as number);
  await expect.poll(() => passing(page), {message: `rows passing the filter (remembered: ${want})`}).toBe(want);
}, {description: 'the filter count equals the remembered one'});

export const passingFewer = Then('fewer rows than remembered should pass the filter', async (page: Page) => {
  const want = await page.evaluate(() => (window as any).__nxPassing as number);
  await expect.poll(() => passing(page), {message: `rows passing the filter (remembered: ${want})`}).toBeLessThan(want);
}, {description: 'the filter count is below the remembered one'});

export const passingMore = Then('more rows than remembered should pass the filter', async (page: Page) => {
  const want = await page.evaluate(() => (window as any).__nxPassing as number);
  await expect.poll(() => passing(page), {message: `rows passing the filter (remembered: ${want})`}).toBeGreaterThan(want);
}, {description: 'the filter count is above the remembered one'});

export const passingNotMore = Then('no more rows than remembered should pass the filter', async (page: Page) => {
  const want = await page.evaluate(() => (window as any).__nxPassing as number);
  await expect.poll(() => passing(page), {message: `rows passing the filter (remembered: ${want})`}).toBeLessThanOrEqual(want);
}, {description: 'the filter count is at most the remembered one'});

// --- the Scaffold Tree filter card (Chem) -----------------------------------------------------------

/* Chem's Scaffold Tree filter hosts a Scaffold Tree viewer that no view owns: its root carries no
   viewer name, so the viewer runtime cannot reach it from the page. Its readings and node areas
   are read through the filter the panel of the view in front holds. */
type TreeStatus = {values: Record<string, unknown>; areas: Record<string, {x: number; y: number; width: number; height: number}>};

async function scaffoldFilterStatus(page: Page): Promise<TreeStatus> {
  return page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    if (!tv?.root.querySelector('[name="viewer-Filters"]'))
      throw new Error('the view in front has no filter panel');
    const filter = (tv.getFiltersGroup().filters as any[]).find((f) => f?.viewer && 'checkedNodes' in f.viewer);
    if (!filter)
      throw new Error('the filter panel holds no Scaffold Tree filter');
    const status = filter.viewer.getWidgetStatus();
    const r = filter.viewer.root.getBoundingClientRect();
    const areas: Record<string, {x: number; y: number; width: number; height: number}> = {};
    for (const [name, a] of Object.entries(status.hitAreas ?? {}) as [string, any][])
      areas[name] = {x: r.left + a.x, y: r.top + a.y, width: a.width, height: a.height};
    return {values: status.values ?? {}, areas};
  });
}

export const scaffoldFilterReading = Then('the {string} reading of the scaffold tree filter should be {int}', async (page: Page, name: string, value: number) => {
  await expect.poll(async () => Number((await scaffoldFilterStatus(page)).values[name] ?? NaN),
    {message: `the "${name}" reading of the scaffold tree filter`}).toBe(value);
}, {description: 'a reading of the Scaffold Tree viewer inside the filter card ("nodes", "checked nodes", "colored nodes", "rows kept")'});

export const scaffoldFilterReadingAtLeast = Then('the {string} reading of the scaffold tree filter should be at least {int}', async (page: Page, name: string, value: number) => {
  await expect.poll(async () => Number((await scaffoldFilterStatus(page)).values[name] ?? NaN),
    {message: `the "${name}" reading of the scaffold tree filter`}).toBeGreaterThanOrEqual(value);
}, {description: 'as above, a lower bound'});

export const clickScaffoldFilterArea = When('user clicks on the {string} area of the scaffold tree filter', async (page: Page, area: string) => {
  let box: {x: number; y: number; width: number; height: number} | undefined;
  await expect.poll(async () => (box = (await scaffoldFilterStatus(page)).areas[area]) !== undefined,
    {timeout: pollMs(5000), message: `the "${area}" area of the scaffold tree filter`}).toBe(true);
  await page.mouse.click(box!.x + box!.width / 2, box!.y + box!.height / 2);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'a click in the middle of a node area the tree reports ("checkbox of node 1")'});

export const scaffoldFilterReadingText = Then('the {string} reading of the scaffold tree filter should be {string}', async (page: Page, name: string, value: string) => {
  await expect.poll(async () => String((await scaffoldFilterStatus(page)).values[name] ?? ''),
    {message: `the "${name}" reading of the scaffold tree filter`}).toBe(value);
}, {description: 'a text reading of the Scaffold Tree viewer inside the filter card ("bit operation")'});

/** The card's toolbar holds the AND / OR choice that combines the checked scaffolds. */
export const setScaffoldBitOperation = When('user sets the scaffold tree filter to combine the checked scaffolds with {string}', async (page: Page, op: string) => {
  const select = page.locator('[name="viewer-Filters"]').filter({visible: true}).first()
    .locator('.d4-filter-element[data-source="Chem:Scaffold Tree Filter"] .chem-scaffold-tree-toolbar select').first();
  await expect(select, 'the AND / OR choice of the scaffold tree filter').toBeVisible({timeout: 5000});
  await select.selectOption({label: op});
  await expect(select).toHaveValue(op, {timeout: 5000});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the AND / OR choice in the toolbar of the Scaffold Tree card, read back'});

/** A grid cell is named by its table row ("cell 962 of Structure"), and which rows a filtered view
 * draws depends on the filter: of the cells drawn in the column, the one whose value is longest is
 * taken — for a molecule, the largest drawn one, which the fewest of the others contain. */
export const pickFromLongestCellMenu = When('user picks {string} from the context menu of the drawn cell of {string} column with the longest value',
  async (page: Page, path: string, column: string) => {
    const grid = el('grid');
    let name = '';
    let box: {x: number; y: number; width: number; height: number} | undefined;
    await expect.poll(async () => {
      const areas = await viewers.hitAreas(page, grid);
      const cells = Object.keys(areas).filter((k) => k.endsWith(` of ${column}`) && /^cell \d+ of /.test(k));
      const lengths = await page.evaluate(([c, rows]) => {
        const df = (window as any).grok.shell.tv.dataFrame;
        return rows.map((r) => String(df.get(c, r - 1) ?? '').length);
      }, [column, cells.map((k) => Number(/^cell (\d+) /.exec(k)![1]))] as [string, number[]]);
      name = cells.map((k, i) => ({k, n: lengths[i]})).sort((a, b) => b.n - a.n)[0]?.k ?? '';
      box = name ? areas[name] : undefined;
      return name !== '';
    }, {timeout: pollMs(5000), message: `a drawn cell of the "${column}" column in the grid`}).toBe(true);
    await page.evaluate(([c, r]) => {
      (window as any).__nxPicked = {column: c, value: (window as any).grok.shell.tv.dataFrame.get(c, r - 1)};
    }, [column, Number(/^cell (\d+) /.exec(name)![1])] as [string, number]);
    await page.mouse.click(box!.x + box!.width / 2, box!.y + box!.height / 2, {button: 'right'});
    await viewers.pickMenuPath(page, path);
  }, {tier: 'ui', description: 'right-clicks, of the cells the grid draws in that column, the one with the longest value, and picks the path in its menu'});

/** A structure card that holds a molecule draws it on a canvas, which opens the sketcher on a click. */
export const clickCardStructure = When('user clicks on the structure drawn in the {string} filter card', async (page: Page, caption: string) => {
  const canvas = page.locator('[name="viewer-Filters"]').filter({visible: true}).first()
    .locator(`[name="filter-card-${caption}"] canvas`).filter({visible: true}).first();
  await expect(canvas, `the structure drawn in the "${caption}" filter card`).toBeVisible({timeout: 5000});
  await canvas.click();
}, {tier: 'ui', description: 'a click on the molecule the substructure card shows'});

/** Chem's Rendering pane of a molecule column sits inside the Chemistry pane and builds its inputs
 * after it opens; a pane header says it is open with the "expanded" class. */
export const setRenderingFilterType = When('user sets Filter type to {string} in the Rendering pane of the context panel', async (page: Page, type: string) => {
  const panel = page.locator('.grok-prop-panel').filter({visible: true}).first();
  for (const section of ['Chemistry', 'Rendering']) {
    const header = panel.locator(`[name="div-section--${section}"]`).first();
    await expect(header, `the ${section} pane of the context panel`).toBeVisible({timeout: 5000});
    if (!/\bexpanded\b/.test(await header.getAttribute('class') ?? ''))
      await header.click();
    await expect(header).toHaveClass(/\bexpanded\b/, {timeout: 5000});
  }
  const select = panel.locator('[name="pane-Rendering"] select[name="input-Filter-type"]').first();
  await expect(select, 'the Filter type input of the Rendering pane').toBeVisible({timeout: pollMs(30000)});
  await select.selectOption({label: type});
  await expect(select).toHaveValue(type, {timeout: 5000});
}, {tier: 'ui', description: 'opens Chemistry > Rendering where it is closed, then the Filter type choice, waited for as long as the pane takes to build (up to 30 s), read back'});

/** What each view's filter panel filters by: per view, the cards that filter and their summaries,
 * read from the panel's own readings — the state a saved project has to bring back. */
async function panelStates(page: Page): Promise<Record<string, string>> {
  return page.evaluate(() => {
    const w = window as any;
    const out: Record<string, string> = {};
    for (const tv of Array.from(w.grok.shell.tableViews as any[])) {
      const panel = Array.from((tv.root as HTMLElement).querySelectorAll('[name="viewer-Filters"]'))[0];
      if (!panel) {
        out[tv.name] = 'no panel';
        continue;
      }
      const v = w.__bdd?.viewerOf(panel) ?? w.DG.Widget.find(panel);
      const values = v?.getWidgetStatus?.()?.values ?? {};
      const cards = String(values['cards'] ?? '').split(', ').filter((c) => values[`filtering of ${c}`] === true);
      // a Scaffold Tree card has no caption in the panel's readings: its tree says what it filters by
      const trees = (tv.getFiltersGroup().filters as any[]).filter((f) => f?.viewer && 'checkedNodes' in f.viewer)
        .map((f) => { const v = f.viewer.getWidgetStatus().values; return `scaffold tree [${v['checked nodes']} of ${v['nodes']} checked, ${v['bit operation']}]`; });
      out[tv.name] = `${values['active'] === false ? 'off: ' : ''}${[...cards.map((c) => `${c} [${values[`summary of ${c}`] ?? ''}]`), ...trees].join('; ')}`;
    }
    const df = w.grok.shell.tv?.dataFrame;
    out['(table filters)'] = df ? Array.from(df.rows.filters as Iterable<string>).map(String).sort().join(' | ') : '';
    return out;
  });
}

export const rememberPanelStates = When('user remembers what the filter panel of every view filters by', async (page: Page) => {
  await viewers.settleAll(page);
  const states = await panelStates(page);
  await page.evaluate((s) => { (window as any).__nxPanels = s; }, states);
}, {tier: 'api', description: 'per table view, the cards that filter and their summaries'});

export const panelStatesAsRemembered = Then('the filter panel of every view should filter by what was remembered', async (page: Page) => {
  const want = await page.evaluate(() => (window as any).__nxPanels);
  await expect.poll(() => panelStates(page), {message: 'what the filter panel of each view filters by'}).toEqual(want);
}, {description: 'every table view has the same filtering cards with the same summaries as remembered'});

// --- the chain's five projects ------------------------------------------------------------------------

/* The listing of every project is slow on a shared stand (minutes on dev): the chain lists once
   before it starts and once to prove the projects gone; in between it finds each project by name. */
const RUN_SUFFIX = /-(\d{13,}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$/;
const namesOf = (list: string): string[] => list.split(/\s*,\s*/).filter(Boolean);

async function deleteProjectsById(page: Page, ids: string[]): Promise<void> {
  await page.evaluate(async (list) => {
    const w = window as any;
    for (const id of list) {
      const project = await w.grok.dapi.projects.filter(`id = "${id}"`).include('children').first();
      if (!project)
        continue;
      for (const child of project.children) {
        const source = child instanceof w.DG.TableInfo ? w.grok.dapi.tables : child instanceof w.DG.ViewInfo ? w.grok.dapi.views : null;
        if (source)
          await source.delete(child);
      }
      await w.grok.dapi.projects.delete(project);
    }
  }, ids);
}

async function projectIdsByName(page: Page, names: string[]): Promise<string[]> {
  return page.evaluate(async (list) => {
    const ids: string[] = [];
    for (const n of list)
      for (const p of await (window as any).grok.dapi.projects.filter(`friendlyName = "${n}"`).list())
        ids.push(String(p.id));
    return ids;
  }, names);
}

/** Every project on the server, by one listing: name, friendly name, id and age. */
async function allProjects(page: Page): Promise<{id: string; name: string; friendlyName: string; age: number}[]> {
  return page.evaluate(async () => (await (window as any).grok.dapi.projects.list() as any[])
    .map((p) => ({id: String(p.id), name: String(p.name), friendlyName: String(p.friendlyName),
      age: Date.now() - (Number(p.createdOn?.valueOf?.() ?? 0) || Date.now())})));
}

export const chainProjectsCleaned = Given('the projects {string} are deleted now and when the feature ends', async (page: Page, list: string) => {
  const names = namesOf(list);
  const families = names.filter((n) => RUN_SUFFIX.test(n)).map((n) => n.replace(RUN_SUFFIX, ''));
  const leftovers = (await allProjects(page)).filter((p) => names.includes(p.friendlyName) ||
    (p.age > 60 * 60 * 1000 && RUN_SUFFIX.test(p.friendlyName) && families.includes(p.friendlyName.replace(RUN_SUFFIX, ''))));
  await deleteProjectsById(page, leftovers.map((p) => p.id));
  atFeatureEnd(page, async () => deleteProjectsById(page, await projectIdsByName(page, names)));
}, {tier: 'api', description: 'one listing of the server\'s projects: those of these names, and those of their families (the same name with another run\'s suffix) older than an hour, are deleted with their tables and views; at feature end each name is looked up and deleted again'});

export const deleteChainProjects = When('user deletes the projects {string}', async (page: Page, list: string) => {
  await deleteProjectsById(page, await projectIdsByName(page, namesOf(list)));
}, {tier: 'api', description: 'each project looked up by its friendly name and deleted with its tables and views'});

export const noChainProjects = Then('none of the projects {string} should be on the server', async (page: Page, list: string) => {
  const names = namesOf(list);
  await expect.poll(async () => (await allProjects(page)).filter((p) => names.includes(p.friendlyName) || names.includes(p.name)).map((p) => p.friendlyName),
    {timeout: pollMs(60000), message: 'projects of these names in the listing of every project'}).toEqual([]);
}, {description: 'read from the listing of every project, which a filter cannot hide an entity from'});

/** The ribbon's Save dialog uploads after OK and says so in the task bar. */
export const saveFromDialog = When('user clicks on OK in the Save project dialog and the project uploads', async (page: Page) => {
  const dialog = page.locator('[name="dialog-Save-project"]').filter({visible: true}).last();
  await gestures.click(page, el('OK button in "Save project" dialog'));
  await expect(dialog, 'the Save project dialog').toBeHidden({timeout: pollMs(60000)});
  await expect(page.locator('.d4-task-bar-entry', {hasText: 'Uploading'}), 'the Uploading entry of the task bar').toHaveCount(0, {timeout: pollMs(120000)});
}, {tier: 'ui', description: 'OK, then the dialog closes and the task bar\'s Uploading entry is gone'});

// --- rows that links carry --------------------------------------------------------------------------

/** The rows of a table whose key tuple is the key tuple of some selected row of another table —
 * what a link keyed that way should leave there — against a bitset of the first table. */
async function linkedRowsCheck(page: Page, table: string, bitset: 'filter' | 'selection', source: string, keys: string, sourceKeys: string): Promise<string> {
  return page.evaluate(([t, b, st, k, sk]) => {
    const tables = (window as any).grok.shell.tables as any[];
    const df = tables.find((x) => x.name === t);
    const src = tables.find((x) => x.name === st);
    if (!df || !src)
      return `no table "${df ? st : t}" is open`;
    const cols = (d: any, names: string): any[] => names.split(/\s*,\s*/).map((n) => {
      const c = d.col(n);
      if (!c)
        throw new Error(`"${d.name}" has no "${n}" column`);
      return c;
    });
    const kc = cols(df, k);
    const skc = cols(src, sk);
    const tuples = new Set<string>();
    for (let i = 0; i < src.rowCount; i++)
      if (src.selection.get(i))
        tuples.add(JSON.stringify(skc.map((c) => c.get(i))));
    const bits = df[b];
    let expected = 0;
    let wrong = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const hit = tuples.has(JSON.stringify(kc.map((c) => c.get(i))));
      if (hit)
        expected++;
      if (hit !== bits.get(i))
        wrong++;
    }
    return wrong === 0 ? `${expected}` : `${bits.trueCount} rows, ${expected} match the ${tuples.size} selected keys of "${st}", ${wrong} rows differ`;
  }, [table, bitset, source, keys, sourceKeys] as [string, 'filter' | 'selection', string, string, string]);
}

async function expectLinked(page: Page, table: string, bitset: 'filter' | 'selection', source: string, keys: string, sourceKeys: string): Promise<void> {
  let got = '';
  await expect.poll(async () => /^\d+$/.test(got = await linkedRowsCheck(page, table, bitset, source, keys, sourceKeys)) && got !== '0',
    {message: `the ${bitset} of "${table}" against the selection of "${source}"`}).toBe(true).catch(() => {
    throw new Error(`the ${bitset} of "${table}": ${got === '0' ? 'no row matches the selection, which proves nothing' : got}`);
  });
}

export const filterMatchesLinkedSelection = Then('the rows of table {string} that pass the filter should be exactly those matching the rows selected in table {string} on {string} = {string}',
  (page: Page, table: string, source: string, keys: string, sourceKeys: string) => expectLinked(page, table, 'filter', source, keys, sourceKeys),
  {description: 'row by row: a row passes exactly when its key columns (comma-separated) hold the keys of some row selected in the other table; at least one row must match'});

export const selectionMatchesLinkedSelection = Then('the rows selected in table {string} should be exactly those matching the rows selected in table {string} on {string} = {string}',
  (page: Page, table: string, source: string, keys: string, sourceKeys: string) => expectLinked(page, table, 'selection', source, keys, sourceKeys),
  {description: 'row by row: a row is selected exactly when its key columns hold the keys of some row selected in the other table; at least one row must match'});

// --- the grid's pinned columns -----------------------------------------------------------------------

export const gridPins = Then('the grid should pin the columns {string}', async (page: Page, list: string) => {
  const want = namesOf(list);
  await expect.poll(() => viewers.onViewer(page, el('grid'), (e) => {
    const w = (window as any).__bdd.viewerOf(e);
    const order = String(w.getWidgetStatus()?.values?.['column order'] ?? '').split(', ');
    return order.slice(0, Number(w.props.frozenColumns) - 1);
  }), {message: 'the columns left of the frozen line (Frozen Columns counts the row header)'}).toEqual(want);
}, {description: 'the first columns of the grid\'s "column order" reading, as many as Frozen Columns keeps (less the row header), in order'});

// --- what a viewer draws of its formula lines --------------------------------------------------------

export const drawsFormulaLines = Then('{widget} should draw {int} formula line(s)', async (page: Page, target: ElementRef, count: number) => {
  let names: string[] = [];
  await expect.poll(async () => (names = Object.keys(await viewers.hitAreas(page, target)).filter((k) => /^formula (line|band) /.test(k))).length,
    {message: `formula line areas of ${target.phrase}`}).toBe(count).catch(() => {
    throw new Error(`${target.phrase} draws ${names.length} formula line(s), not ${count}: ${names.join(', ') || 'none'}`);
  });
}, {description: 'the "formula line <title>" / "formula band <title>" hit areas, which a viewer reports only for an item it drew on screen'});

export const propertyLacks = Then('{string} property of {widget} should not contain {string}', async (page: Page, caption: string, target: ElementRef, text: string) => {
  const value = await viewers.readProperty(page, target, caption);
  expect(value, `"${caption}" of ${target.phrase}`).not.toBe('');
  expect(value, `"${caption}" of ${target.phrase}`).not.toContain(text);
}, {description: 'the property is set and its value as text holds no such substring'});

export const tableTagLacks = Then('the {string} tag of the table should not contain {string}', async (page: Page, tag: string, text: string) => {
  const value = await page.evaluate((t) => String((window as any).grok.shell.tv?.dataFrame?.getTag(t) ?? ''), tag);
  expect(value, `the "${tag}" tag of the current table`).not.toBe('');
  expect(value, `the "${tag}" tag of the current table`).not.toContain(text);
}, {description: 'the tag of the table in front is set and holds no such substring'});

export const formulaLineRanges = Then('the formula lines of {widget} should hold {string} over the ranges {string}', async (page: Page, target: ElementRef, formula: string, ranges: string) => {
  const want = ranges.split(/\s*;\s*/).filter(Boolean).sort();
  await expect.poll(async () => {
    const raw = await viewers.readProperty(page, target, 'formulaLines');
    return (JSON.parse(raw || '[]') as any[]).filter((l) => l.formula === formula).map((l) => `${l.min ?? ''}..${l.max ?? ''}`).sort();
  }, {message: `the items of "formulaLines" of ${target.phrase} with that formula, as min..max`}).toEqual(want);
}, {description: 'the items of the viewer\'s formulaLines look with exactly this formula, their min..max pairs (";"-separated, any order)'});

// --- a remembered count by name ---------------------------------------------------------------------

export const rememberPassingAs = When('user remembers how many rows pass the filter as {string}', async (page: Page, name: string) => {
  await viewers.settleAll(page);
  await page.evaluate((n) => {
    const w = window as any;
    w.__nxCounts = {...(w.__nxCounts ?? {}), [n]: {view: String(w.grok.shell.tv.name), count: w.grok.shell.tv.dataFrame.filter.trueCount}};
  }, name);
}, {tier: 'api', description: 'the filter count of the table in front, kept under a name with the view it was read in'});

export const passingAsRememberedAs = Then('as many rows as remembered as {string} should pass the filter', async (page: Page, name: string) => {
  const want = await page.evaluate((n) => ((window as any).__nxCounts ?? {})[n], name) as {view: string; count: number} | undefined;
  if (!want)
    throw new Error(`no count was remembered as "${name}"`);
  await expect.poll(() => page.evaluate((v) => {
    const tv = Array.from((window as any).grok.shell.tableViews as any[]).find((x) => x.name === v);
    return tv ? tv.dataFrame.filter.trueCount : `no "${v}" view is open`;
  }, want.view), {message: `rows passing the filter of the "${want.view}" view's table (remembered as "${name}": ${want.count})`}).toBe(want.count);
}, {description: 'the filter count of the table of the view the count was read in equals the one remembered under that name'});

// --- Chem filters that compute ------------------------------------------------------------------------

/** A substructure card searches in the background ("searching of <col>"), and a Scaffold Tree counts
 * the hits of its nodes after it loads (-1 until then): a count read before both are done is the
 * count of a filter not applied yet. */
export const chemFiltersDone = Then('the Chem filters of every view should have finished computing', async (page: Page) => {
  let pending: string[] = [];
  await expect.poll(async () => (pending = await page.evaluate(() => {
    const w = window as any;
    const out: string[] = [];
    for (const tv of Array.from(w.grok.shell.tableViews as any[])) {
      const panel = (tv.root as HTMLElement).querySelector('[name="viewer-Filters"]');
      if (!panel)
        continue;
      const values = w.__bdd?.viewerOf(panel)?.getWidgetStatus?.()?.values ?? {};
      for (const [k, v] of Object.entries(values))
        if (k.startsWith('searching of ') && v === true)
          out.push(`${tv.name}: ${k}`);
      for (const f of tv.getFiltersGroup().filters as any[]) {
        if (!(f?.viewer && 'checkedNodes' in f.viewer))
          continue;
        const tree = f.viewer.getWidgetStatus().values;
        for (let i = 1; i <= Number(tree['nodes'] ?? 0); i++)
          if (tree[`checked of node ${i}`] === true && Number(tree[`hits of node ${i}`]) < 0)
            out.push(`${tv.name}: hits of checked node ${i}`);
      }
    }
    return out;
  })).length, {timeout: pollMs(60000), message: 'Chem filters still computing'}).toBe(0).catch(() => {
    throw new Error(`Chem filters still computing: ${pending.join(', ')}`);
  });
  await viewers.settleAll(page);
}, {description: 'no substructure card of any view is searching and every checked Scaffold Tree node has counted its hits (up to a minute)'});

export const scaffoldTreeUncolored = Then('no node of the scaffold tree filter should be colored', async (page: Page) => {
  let values: Record<string, unknown> = {};
  await expect.poll(async () => {
    values = (await scaffoldFilterStatus(page)).values;
    const n = Number(values['nodes'] ?? 0);
    return n > 0 && Array.from({length: n}, (_, i) => Number(values[`hits of node ${i + 1}`])).every((h) => h >= 0);
  }, {timeout: pollMs(60000), message: 'every node of the scaffold tree filter has counted its hits'}).toBe(true);
  const n = Number(values['nodes']);
  const colored = Array.from({length: n}, (_, i) => `node ${i + 1}: ${values[`color of node ${i + 1}`] ?? ''}`).filter((s) => !s.endsWith(': '));
  expect(colored, `colored nodes of the scaffold tree filter (its "colored nodes" reading: ${values['colored nodes']})`).toEqual([]);
  expect(Number(values['colored nodes']), 'the "colored nodes" reading').toBe(0);
}, {description: 'once every node has counted its hits (the tree is built), each node\'s color reading is empty and the "colored nodes" reading is 0'});

// --- a calculated column's formula ---------------------------------------------------------------------

export const formulaColumnsExist = Then('every column the formula of {string} column refers to should exist', async (page: Page, column: string) => {
  const missing = await page.evaluate((c) => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    const formula = String(df?.col(c)?.getTag('formula') ?? '');
    if (!formula)
      return [`"${c}" has no formula`];
    return Array.from(formula.matchAll(/\$\{([^}]+)\}/g), (m) => m[1]).filter((n) => !df.col(n));
  }, column);
  expect(missing, `columns the formula of "${column}" names that the table does not have`).toEqual([]);
}, {description: 'every ${name} of the column\'s formula tag is a column of the table'});

/** The Save project dialog draws its preview of the view as it opens: the window for that preview's
 * console noise is opened with the dialog. */
export const openSaveDialog = When('user opens the Save project dialog from the ribbon', async (page: Page) => {
  projectSaveWindow(page, true);
  await gestures.click(page, el('Save button'));
  await expect(page.locator('[name="dialog-Save-project"]').filter({visible: true}), 'the Save project dialog').toHaveCount(1, {timeout: pollMs(15000)});
}, {tier: 'ui', description: 'the ribbon\'s Save button; done when the Save project dialog is shown. The preview\'s "cloned iframe" console message is ignored from here until the next project opens'});

/** A formula line keeps the title it was given when it was made; what it computes is its formula. */
const formulasOf = (raw: string): string[] => (JSON.parse(raw || '[]') as any[]).map((l) => String(l.formula ?? ''));

export const viewerFormulasLack = Then('no formula line of {widget} should name {string}', async (page: Page, target: ElementRef, text: string) => {
  const formulas = formulasOf(await viewers.readProperty(page, target, 'formulaLines'));
  expect(formulas.length, `formula lines of ${target.phrase}`).toBeGreaterThan(0);
  expect(formulas.filter((f) => f.includes(text)), `formulas of ${target.phrase} that name "${text}"`).toEqual([]);
}, {description: 'the formulas of the viewer\'s formula lines (not their titles), of which there must be some'});

export const tableFormulasLack = Then('no formula line of the table should name {string}', async (page: Page, text: string) => {
  const formulas = formulasOf(await page.evaluate(() => String((window as any).grok.shell.tv?.dataFrame?.getTag('.formula-lines') ?? '')));
  expect(formulas.length, 'formula lines of the current table').toBeGreaterThan(0);
  expect(formulas.filter((f) => f.includes(text)), `formulas of the table that name "${text}"`).toEqual([]);
}, {description: 'the formulas of the table\'s ".formula-lines" (not their titles), of which there must be some'});

export const rememberScaffoldReading = When('user remembers the {string} reading of the scaffold tree filter', async (page: Page, name: string) => {
  await viewers.settleAll(page);
  const value = Number((await scaffoldFilterStatus(page)).values[name]);
  await page.evaluate(([n, v]) => { (window as any).__nxTree = {...((window as any).__nxTree ?? {}), [n]: v}; }, [name, value] as [string, number]);
}, {tier: 'api', description: 'a numeric reading of the Scaffold Tree filter, for a comparison after a change'});

export const scaffoldReadingLower = Then('the {string} reading of the scaffold tree filter should be lower than remembered', async (page: Page, name: string) => {
  const want = await page.evaluate((n) => ((window as any).__nxTree ?? {})[n], name);
  if (typeof want !== 'number')
    throw new Error(`the "${name}" reading of the scaffold tree filter was not remembered`);
  await expect.poll(async () => Number((await scaffoldFilterStatus(page)).values[name]),
    {message: `the "${name}" reading of the scaffold tree filter (remembered: ${want})`}).toBeLessThan(want);
}, {description: 'the reading is below the one remembered'});

// --- viewers stacked as tabs of one panel --------------------------------------------------------------

/* Dragging a viewer by its title bar over another shows dock-spawn's compass over that viewer; its
   middle item (`.dock-wheel-fill`) puts the dragged viewer into the other's panel as a tab. A tabbed
   panel is a `.dock-container-fill` with one `view-handle: <title>` tab handle per viewer, and only
   the selected tab's viewer is in the page. */
export const stackViewer = When('user stacks {widget} onto {widget} as a tab', async (page: Page, target: ElementRef, other: ElementRef) => {
  const panel = (await viewers.viewerLocator(page, target)).locator('xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]');
  const title = panel.locator('.panel-titlebar').first();
  const t = await title.boundingBox();
  if (!t)
    throw new Error(`${target.phrase} has no title bar to drag`);
  const start = {x: t.x + Math.min(30, t.width / 4), y: t.y + t.height / 2};
  await page.mouse.move(start.x, start.y);
  await page.mouse.down();
  await page.mouse.move(start.x + 20, start.y + 20, {steps: 4});
  const fill = page.locator('.dock-wheel-base .dock-wheel-item.dock-wheel-fill').filter({visible: true}).first();
  try {
    // the dragged panel leaves the layout once the drag starts, and the others take its room:
    // the target is measured after that
    const over = await (await viewers.viewerLocator(page, other)).boundingBox();
    if (!over)
      throw new Error(`${other.phrase} has no box to drop on`);
    await page.mouse.move(over.x + over.width / 2, over.y + over.height / 2, {steps: 8});
    await fill.waitFor({timeout: pollMs(5000)}).catch(() => {
      throw new Error(`dragging ${target.phrase} over ${other.phrase} showed no dock compass`);
    });
    const box = (await fill.boundingBox())!;
    const inside = box.x + box.width / 2 > over.x && box.x + box.width / 2 < over.x + over.width &&
      box.y + box.height / 2 > over.y && box.y + box.height / 2 < over.y + over.height;
    if (!inside)
      throw new Error(`the dock compass did not open over ${other.phrase}`);
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2, {steps: 6});
    await expect(fill, 'the middle item of the dock compass under the pointer').toHaveClass(/dock-wheel-fill-icon-hover/, {timeout: pollMs(3000)});
  }
  catch (e) {
    await page.mouse.up();
    throw e;
  }
  await page.mouse.up();
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'drags the viewer by its title bar onto the middle item of the dock compass over the other viewer, which makes it a tab of that viewer\'s panel'});

/** The tabs of the tabbed panels of the view in front, by panel. */
async function tabPanels(page: Page): Promise<{tabs: string[]; selected: string}[]> {
  return page.evaluate(() => {
    const root = (window as any).grok.shell.v?.root as HTMLElement | undefined;
    return Array.from(root?.querySelectorAll('.dock-container-fill') ?? []).map((c) => {
      const handles = Array.from(c.querySelectorAll('.tab-handle[name^="view-handle: "]'))
        .filter((h) => h.closest('.dock-container-fill') === c);
      return {tabs: handles.map((h) => h.getAttribute('name')!.slice('view-handle: '.length)),
        selected: handles.filter((h) => h.classList.contains('tab-handle-selected')).map((h) => h.getAttribute('name')!.slice('view-handle: '.length))[0] ?? ''};
    }).filter((p) => p.tabs.length > 1);
  });
}

export const viewersTabbed = Then('the viewers {string} should be the tabs of one panel', async (page: Page, list: string) => {
  const want = namesOf(list).sort();
  await expect.poll(async () => (await tabPanels(page)).map((p) => [...p.tabs].sort().join(', ')),
    {message: 'the tabbed panels of the current view, by their tabs'}).toContain(want.join(', '));
}, {description: 'one tabbed panel of the view in front holds exactly these viewers (their titles, comma-separated) as its tabs'});

export const noTabbedPanel = Then('no panel of the current view should hold tabs', async (page: Page) => {
  await expect.poll(async () => (await tabPanels(page)).map((p) => p.tabs.join(', ')), {message: 'the tabbed panels of the current view'}).toEqual([]);
}, {description: 'every viewer of the view in front has a panel of its own'});

export const switchTab = When('user switches the tabbed panel to the {string} tab', async (page: Page, title: string) => {
  const handle = page.locator(`.dock-container-fill .tab-handle[name="view-handle: ${title}"]`).filter({visible: true}).first();
  await expect(handle, `the "${title}" tab of a tabbed panel`).toBeVisible({timeout: 5000});
  await handle.click();
  await expect.poll(async () => (await tabPanels(page)).some((p) => p.selected === title), {message: `the "${title}" tab selected`}).toBe(true);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'a click on the tab handle of that title; done when it is the selected tab and the viewers have settled'});


export const rowsContainPicked = Then('every row that passes the filter should contain the molecule of the cell picked', async (page: Page) => {
  const r = await page.evaluate(async () => {
    const w = window as any;
    const picked = w.__nxPicked;
    if (!picked)
      return {error: 'no cell was picked'};
    const df = w.grok.shell.tv.dataFrame;
    const rdkit = await w.grok.functions.call('Chem:getRdKitModule');
    const query = rdkit.get_qmol(String(picked.value));
    const rows = Array.from(df.filter.getSelectedIndexes() as Iterable<number>);
    const without: number[] = [];
    try {
      for (const i of rows) {
        const mol = rdkit.get_mol(String(df.get(picked.column, i)));
        if (mol.get_substruct_match(query) === '{}')
          without.push(i + 1);
        mol.delete();
      }
    }
    finally {
      query.delete();
    }
    return {passing: rows.length, without};
  });
  if ('error' in r)
    throw new Error(r.error);
  expect(r.passing, 'rows passing the filter').toBeGreaterThan(0);
  expect(r.without, `rows passing the filter (of ${r.passing}) whose molecule does not contain the picked one`).toEqual([]);
}, {description: 'RDKit substructure match: every row that passes holds a molecule that contains the molecule of the cell the context menu was opened on'});
