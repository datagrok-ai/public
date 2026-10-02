/* The steps only the NX chain (features/viewers/nx) needs: a viewer's rows against its own formula
   filter, Chem's Scaffold Tree filter card and substructure card, and the molecule a cell's menu took
   as a filter. Everything else the chain does — Link Tables, the Save project dialog, view tabs,
   layouts kept by name, formula lines, tabbed viewers, filter counts across views — is library
   vocabulary. */
import {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {cssString, el, type ElementRef, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

/** The rows of a table that pass its filter (what a link brought) and a viewer formula filter of the
 * "<text column> is <value> and <numeric column> is below <n>" shape — what a viewer bound to that
 * table with that formula should draw. */
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

// --- the Scaffold Tree filter card (Chem) -----------------------------------------------------------

/* Chem's Scaffold Tree filter hosts a Scaffold Tree viewer that no view owns: its root is a plain box
   with no viewer name, so the viewer runtime (and every {widget} step) cannot reach it from the page.
   Its readings and node areas are read through the filter the panel of the view in front holds. */
type TreeStatus = {values: Record<string, unknown>; areas: Record<string, {x: number; y: number; width: number; height: number}>;
  missing?: string};

/** The tree's readings and node areas, or why there are none yet: a card still being built is not a
 * failure a poll should stop on. */
async function scaffoldFilterStatus(page: Page): Promise<TreeStatus> {
  return page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    if (!tv?.root.querySelector('[name="viewer-Filters"]'))
      return {values: {}, areas: {}, missing: 'the view in front has no filter panel'};
    const filter = (tv.getFiltersGroup().filters as any[]).find((f) => f?.viewer && 'checkedNodes' in f.viewer);
    if (!filter)
      return {values: {}, areas: {}, missing: 'the filter panel holds no Scaffold Tree filter'};
    const status = filter.viewer.getWidgetStatus();
    const r = filter.viewer.root.getBoundingClientRect();
    const areas: Record<string, {x: number; y: number; width: number; height: number}> = {};
    for (const [name, a] of Object.entries(status.hitAreas ?? {}) as [string, any][])
      areas[name] = {x: r.left + a.x, y: r.top + a.y, width: a.width, height: a.height};
    return {values: status.values ?? {}, areas};
  });
}

/** A reading as the poll sees it: the number, or the reason there is none. */
async function treeReading(page: Page, name: string): Promise<number | string> {
  const st = await scaffoldFilterStatus(page);
  return st.missing ?? (st.values[name] === undefined ? `no "${name}" reading` : Number(st.values[name]));
}

export const scaffoldFilterReading = Then('the {string} reading of the scaffold tree filter should be {int}', async (page: Page, name: string, value: number) => {
  await expect.poll(() => treeReading(page, name), {message: `the "${name}" reading of the scaffold tree filter`}).toBe(value);
}, {description: 'a reading of the Scaffold Tree viewer inside the filter card ("nodes", "checked nodes", "colored nodes", "rows kept")'});

export const scaffoldFilterReadingAtLeast = Then('the {string} reading of the scaffold tree filter should be at least {int}', async (page: Page, name: string, value: number) => {
  let seen: number | string = '';
  await expect.poll(async () => typeof (seen = await treeReading(page, name)) === 'number' && seen >= value,
    {message: `the "${name}" reading of the scaffold tree filter, at least ${value}`}).toBe(true).catch(() => {
    throw new Error(`the "${name}" reading of the scaffold tree filter is ${seen}, not at least ${value}`);
  });
}, {description: 'as above, a lower bound'});

export const scaffoldFilterReadingText = Then('the {string} reading of the scaffold tree filter should be {string}', async (page: Page, name: string, value: string) => {
  await expect.poll(async () => {
    const st = await scaffoldFilterStatus(page);
    return st.missing ?? String(st.values[name] ?? '');
  }, {message: `the "${name}" reading of the scaffold tree filter`}).toBe(value);
}, {description: 'a text reading of the Scaffold Tree viewer inside the filter card ("bit operation")'});

export const clickScaffoldFilterArea = When('user clicks on the {string} area of the scaffold tree filter', async (page: Page, area: string) => {
  let box: {x: number; y: number; width: number; height: number} | undefined;
  await expect.poll(async () => (box = (await scaffoldFilterStatus(page)).areas[area]) !== undefined,
    {timeout: pollMs(5000), message: `the "${area}" area of the scaffold tree filter`}).toBe(true);
  await page.mouse.click(box!.x + box!.width / 2, box!.y + box!.height / 2);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'a click in the middle of a node area the tree reports ("checkbox of node 1")'});

/** The card's toolbar holds the AND / OR choice that combines the checked scaffolds. */
export const setScaffoldBitOperation = When('user sets the scaffold tree filter to combine the checked scaffolds with {string}', async (page: Page, op: string) => {
  const select = page.locator('[name="viewer-Filters"]').filter({visible: true}).first()
    .locator('.d4-filter-element[data-source="Chem:Scaffold Tree Filter"] .chem-scaffold-tree-toolbar select').first();
  await expect(select, 'the AND / OR choice of the scaffold tree filter').toBeVisible({timeout: pollMs(5000)});
  await select.selectOption({label: op});
  await expect(select).toHaveValue(op, {timeout: pollMs(5000)});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the AND / OR choice in the toolbar of the Scaffold Tree card, read back'});

const rememberedTreeReadings = new WeakMap<Page, Map<string, number>>();

export const rememberScaffoldReading = When('user remembers the {string} reading of the scaffold tree filter', async (page: Page, name: string) => {
  await viewers.settleAll(page);
  const reading = await treeReading(page, name);
  if (typeof reading !== 'number')
    throw new Error(`the "${name}" reading of the scaffold tree filter cannot be remembered: ${reading}`);
  if (!rememberedTreeReadings.has(page))
    rememberedTreeReadings.set(page, new Map());
  rememberedTreeReadings.get(page)!.set(name, reading);
}, {tier: 'api', description: 'a numeric reading of the Scaffold Tree filter, for a comparison after a change'});

export const scaffoldReadingLower = Then('the {string} reading of the scaffold tree filter should be lower than remembered', async (page: Page, name: string) => {
  const want = rememberedTreeReadings.get(page)?.get(name);
  if (want === undefined)
    throw new Error(`the "${name}" reading of the scaffold tree filter was not remembered`);
  let seen: number | string = '';
  await expect.poll(async () => typeof (seen = await treeReading(page, name)) === 'number' && seen < want,
    {message: `the "${name}" reading of the scaffold tree filter (remembered: ${want})`}).toBe(true).catch(() => {
    throw new Error(`the "${name}" reading of the scaffold tree filter is ${seen}, not lower than the remembered ${want}`);
  });
}, {description: 'the reading is below the one remembered'});

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

// --- a molecule used as a filter (Chem) -----------------------------------------------------------------

/** A structure card that holds a molecule draws it on a canvas, which opens the sketcher on a click. */
export const clickCardStructure = When('user clicks on the structure drawn in the {string} filter card', async (page: Page, caption: string) => {
  const canvas = page.locator('[name="viewer-Filters"]').filter({visible: true}).first()
    .locator(`[name="filter-card-${cssString(caption)}"] canvas`).filter({visible: true}).first();
  await expect(canvas, `the structure drawn in the "${caption}" filter card`).toBeVisible({timeout: pollMs(5000)});
  await canvas.click();
}, {tier: 'ui', description: 'a click on the molecule the substructure card shows'});

const pickedCells = new WeakMap<Page, {column: string; value: string}>();

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
    const value = await page.evaluate(([c, r]) => String((window as any).grok.shell.tv.dataFrame.get(c, r - 1)),
      [column, Number(/^cell (\d+) /.exec(name)![1])] as [string, number]);
    pickedCells.set(page, {column, value});
    await page.mouse.click(box!.x + box!.width / 2, box!.y + box!.height / 2, {button: 'right'});
    await viewers.pickMenuPath(page, path);
  }, {tier: 'ui', description: 'right-clicks, of the cells the grid draws in that column, the one with the longest value, and picks the path in its menu'});

export const rowsContainPicked = Then('every row that passes the filter should contain the molecule of the cell picked', async (page: Page) => {
  const picked = pickedCells.get(page);
  if (!picked)
    throw new Error('no cell was picked — "user picks … from the context menu of the drawn cell …" comes first');
  let r: {passing: number; without: number[]} = {passing: 0, without: []};
  // read until the substructure search has settled; a read that throws (RDKit not loaded) says so
  await expect.poll(async () => {
    r = await rowsWithoutPicked(page, picked);
    return r.passing > 0 && r.without.length === 0;
  }, {message: 'the rows passing the filter against the molecule picked'}).toBe(true).catch((e) => {
    if (!/rows passing the filter against/.test(String(e?.message ?? e)))
      throw e;
    throw new Error(r.passing === 0 ? 'no row passes the filter, which proves nothing' :
      `rows passing the filter (of ${r.passing}) whose molecule does not contain the picked one: ${r.without.slice(0, 20).join(', ')}`);
  });
}, {description: 'RDKit substructure match, polled: every row that passes holds a molecule that contains the molecule of the cell the context menu was opened on, and some row passes'});

async function rowsWithoutPicked(page: Page, picked: {column: string; value: string}): Promise<{passing: number; without: number[]}> {
  return page.evaluate(async (p) => {
    const w = window as any;
    const df = w.grok.shell.tv.dataFrame;
    const rdkit = await w.grok.functions.call('Chem:getRdKitModule');
    const query = rdkit.get_qmol(p.value);
    const rows = Array.from(df.filter.getSelectedIndexes() as Iterable<number>);
    const without: number[] = [];
    try {
      for (const i of rows) {
        const mol = rdkit.get_mol(String(df.get(p.column, i)));
        if (mol.get_substruct_match(query) === '{}')
          without.push(i + 1);
        mol.delete();
      }
    }
    finally {
      query.delete();
    }
    return {passing: rows.length, without};
  }, picked);
}
