/* V-1/V-2 — platform viewers as u2 controls (VIEWERS_AND_WIDGETS.md, viewers/plan.md WO-V4): the
   master–detail sample with a real DG.Grid, FilterGroup and Form over a `sample`-policy source.
   Design time shows the sample rows without a run; in Run mode the frame is the bus — a grid click
   moves a bound field, a filter click narrows the grid and the form together; the panel groups a
   viewer's look by its own categories and a checkbox edit is one patch; a look prop binds two-way
   through the picker; a view close kills every viewer the spec built. */
import {artifact, ok, openApp, shot} from '../local.mjs';
import {APP, bindField, bindThroughPicker, dropControl, dumpViaDialog, expandTree, fieldValue, focusTree,
  gridCellCenter, openSection, openSpec, panel, panelField, pickLeaf, reopenApp, row, selectRow, setField, specErrors,
  spyRuns, statusText, surfaceNode, toMode, viewerAt, waitCaption, waitStatus} from '../lib.mjs';

/** `VIEWER_SAMPLES[0]` (src/dg/viewers/samples.ts), verbatim. */
const SAMPLE = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-func-source', name: 'orders',
    props: {func: 'demoOrders', params: {days: 30}, designData: 'sample', debounce: 50,
      sample: [{orderId: 1001, customer: 'Aspirin Labs', city: 'Kyiv', total: 1240},
        {orderId: 1002, customer: 'Bayer', city: 'Lviv', total: 380},
        {orderId: 1003, customer: 'Roche', city: 'Basel', total: 2150}]}}],
  root: {tag: 'u2-splitter', name: 'masterDetail', props: {sizes: [0.25, 0.45, 0.3]}, children: [
    {tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'},
      props: {columnNames: ['city', 'total']}},
    {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.orders'}, props: {allowEdit: false}},
    {tag: 'u2-panel', name: 'detailPane', children: [
      {tag: 'u2-viewer-form', name: 'form', bind: {table: '$.orders'}},
      {tag: 'u2-form', name: 'detailForm', children: [
        {tag: 'u2-text-input', name: 'customerField', props: {label: 'Customer'},
          bind: {value: '$.orders.currentRow.customer'}},
        {tag: 'u2-number-input', name: 'rowIdx', props: {label: 'Row', mode: 'int'},
          bind: {value: '$.orders.currentRowIdx'}}]}]}]},
}, null, 2);

/** Armed before the spec opens: the DD9 proof is that `demoOrders` never ran at design time. */
export async function fixture(page) {
  await page.locator('.u2-palette input').first().fill('');
  await spyRuns(page);
  await openSpec(page, SAMPLE);
  await waitStatus(page, 'masterDetail');
  await page.waitForTimeout(1500);
  await expandTree(page);
}

const dump = async (page) => JSON.parse(await dumpViaDialog(page));

const findNode = (node, name) => node.name === name ? node :
  (node.children ?? []).map((child) => findNode(child, name)).find((found) => found !== undefined);

/** Design mode puts a glass pane over the canvas, and Playwright refuses `locator.click` on what it
 * covers — the selection is made the way a user makes it, with the pointer; the one-behind context
 * panel is re-asserted like `selectRow` does. */
async function selectOnCanvas(page, name) {
  const click = async () => {
    const box = await surfaceNode(page, name).boundingBox();
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  };
  await click();
  await waitStatus(page, name);
  if (!(await waitCaption(page, name, 1500)).startsWith(`${name} (`)) {
    await click();
    await waitCaption(page, name);
  }
  return statusText(page);
}

/** The categorical filter is itself a grid (canvas): the centre of the name cell that reads `value`,
 * in whichever pane of `[data-u2-name=name]` shows it. */
const filterCategoryCenter = (page, name, value) => page.evaluate(({name, value}) => {
  const root = document.querySelector(`.u2-designer-surface [data-u2-name="${name}"]`);
  for (const el of root.querySelectorAll('.d4-filter-categorical')) {
    const grid = DG.Widget.find(el);
    if (!grid)
      continue;
    for (const column of grid.dataFrame.columns.names()) {
      for (let i = 0; i < grid.dataFrame.rowCount; i++) {
        const cell = grid.cell(column, i);
        if (cell.cell.value !== value)
          continue;
        const b = cell.bounds;
        const r = grid.root.getBoundingClientRect();
        return {x: r.left + b.x + b.width / 2, y: r.top + b.y + b.height / 2, column, row: i};
      }
    }
  }
  return null;
}, {name, value});

const filterPanes = (page) => page.evaluate(() =>
  document.querySelectorAll('.u2-designer-surface [data-u2-name="filters"] .d4-filter-element').length);

async function checkDesignTime(page) {
  const pane = await page.evaluate(() => ({
    titles: [...document.querySelectorAll('.u2-palette .u2-accordion-title')].map((t) => t.textContent),
    viewers: document.querySelectorAll('.u2-palette-item[data-u2-tag^="u2-viewer-"]').length,
    // package-internal descriptors (`_makeInspectorPanel`) are no viewers
    descriptors: DG.WidgetDescriptor.getDescriptors().filter((d) => !d.name.startsWith('_')).length,
  }));
  ok('viewers/1/the-palette-lists-every-platform-viewer-under-viewers',
    pane.titles.includes('Viewers') && pane.viewers > 0 && pane.viewers === pane.descriptors, JSON.stringify(pane));

  const grid = await viewerAt(page, 'grid');
  const ran = await page.evaluate(() => window.__u2Runs);
  await shot(page, 'viewers-1-design');
  ok('viewers/1a/design-time-shows-the-sample-rows-in-a-real-grid-without-running-the-source',
    grid?.type === 'Grid' && grid.rows === 3 && grid.canvases > 0 && await specErrors(page) === 0 &&
    !ran.includes('demoOrders'),
    `grid=${JSON.stringify(grid)} broken=${await specErrors(page)} ran=${JSON.stringify(ran)}`);

  const status = await selectOnCanvas(page, 'grid');
  const shown = await panel(page);
  const allowEdit = Object.values(shown).find((section) => 'allowEdit' in section)?.allowEdit;
  await shot(page, 'viewers-1-grid-selected');
  ok('viewers/1b/a-canvas-click-selects-the-grid-and-the-panel-groups-its-look-by-category',
    /grid · \d+ nodes/.test(status) && ['Data', 'Columns', 'Rows', 'Style'].every((s) => s in shown) &&
    allowEdit === 'false',
    `status="${status}" sections=${JSON.stringify(Object.keys(shown))} allowEdit=${allowEdit}`);
}

async function checkRunMode(page) {
  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  const run = await viewerAt(page, 'grid');
  const at = await gridCellCenter(page, 'grid', 'customer', 2);
  await page.mouse.click(at.x, at.y);
  await page.waitForTimeout(600);
  const clicked = await fieldValue(page, 'customerField');
  await setField(page, 'rowIdx', '0');
  await page.waitForTimeout(400);
  const moved = await viewerAt(page, 'grid');
  const stepped = await fieldValue(page, 'customerField');
  await shot(page, 'viewers-1-run');
  ok('viewers/1c/run-mode-runs-the-source-and-the-frame-carries-the-current-row-both-ways',
    run?.rows === 4 && clicked === 'Roche' && moved?.current === 0 && stepped === 'Aspirin Labs',
    `rows=${run?.rows} click row 2 → "${clicked}"; rowIdx 0 → current=${moved?.current} "${stepped}"`);

  const filters = await viewerAt(page, 'filters');
  const lviv = await filterCategoryCenter(page, 'filters', 'Lviv');
  if (lviv !== null) {
    await page.mouse.click(lviv.x, lviv.y);
    await page.waitForTimeout(800);
  }
  const grid = await viewerAt(page, 'grid');
  const form = await viewerAt(page, 'form');
  await shot(page, 'viewers-1-filtered');
  ok('viewers/1d/a-filter-click-narrows-the-grid-and-the-form-through-the-frame',
    filters?.type === 'Filters' && lviv !== null && grid?.filtered === 1 && form?.filtered === 1,
    `filters=${JSON.stringify(filters)} lviv=${JSON.stringify(lviv)} grid.filtered=${grid?.filtered} ` +
    `form.filtered=${form?.filtered}`);
  await toMode(page, 'Design');
  await page.waitForTimeout(500);
}

/** VP-11 — the one verb on `u2-viewer-filters`. The sample keeps its `columnNames` panes — the
 * platform turns them into live `filters` at build — so the picked column's pane is appended after
 * them (P1), and the written prop carries the existing panes first. */
async function checkAddFilter(page) {
  const panesBefore = await filterPanes(page);
  await row(page, 'filters').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  const verb = page.locator('.u2-menu-item').filter({hasText: /^Add filter for column…$/});
  const offered = await verb.count() === 1;
  if (offered)
    await verb.first().click();
  else
    await page.keyboard.press('Escape');
  const dialog = page.locator('.u2-dialog').filter({hasText: 'Column'});
  const opened = await dialog.first().waitFor({state: 'visible', timeout: 3000}).then(() => true, () => false);
  if (opened) {
    await dialog.first().locator('select').first().selectOption('orderId');
    await dialog.first().locator('button').filter({hasText: /^OK$/i}).first().click();
    await page.waitForTimeout(800);
  }
  const added = (await dump(page)).root.children[0].props?.filters;
  const panesAfter = await filterPanes(page);
  await shot(page, 'viewers-1-add-filter');
  await focusTree(page, 'filters');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(600);
  const undone = (await dump(page)).root.children[0].props?.filters;
  const last = added?.[added.length - 1];
  ok('viewers/1e/add-filter-for-column-writes-the-filters-prop-and-rebuilds-the-panes',
    offered && opened && added?.length === panesBefore + 1 && last?.type === 'histogram' && last?.column === 'orderId' &&
    panesAfter === panesBefore + 1 && undone === undefined && await filterPanes(page) === panesBefore,
    `offered=${offered} dialog=${opened} filters=${JSON.stringify(added)} panes ${panesBefore}→${panesAfter}` +
    `→${await filterPanes(page)} after undo=${JSON.stringify(undone)}`);
}

async function checkPanelEdit(page) {
  await selectRow(page, 'grid');
  await openSection(page, 'Data');
  const box = panelField(page, 'allowEdit');
  await box.check();
  await page.waitForTimeout(800);
  const live = await viewerAt(page, 'grid', 'allowEdit');
  const written = (await dump(page)).root.children[1].props?.allowEdit;
  await shot(page, 'viewers-1-allow-edit');
  ok('viewers/1f/a-checkbox-edit-in-the-panel-reaches-the-viewer-and-the-document',
    live?.props?.allowEdit === true && written === true,
    `viewer.props.allowEdit=${live?.props?.allowEdit} dump.allowEdit=${written}`);
}

/** A look prop as a bind target: a choice input over the frame's columns drives the scatter plot's
 * x axis through `$.plot.xColumnName`, picked — not typed — in the binding picker. */
async function checkBoundLookProp(page) {
  await dropControl(page, 'u2-viewer-scatter-plot', 'scatter', 'detailPane');
  await waitStatus(page, 'viewerScatterPlot1');
  await row(page, 'viewerScatterPlot1').locator('.u2-tree-label').dblclick();
  const editor = page.locator('.u2-tree-rename');
  await editor.fill('plot');
  await editor.press('Enter');
  await page.waitForTimeout(400);
  await selectRow(page, 'plot');
  await openSection(page, 'Bindings');
  // a viewer's Bindings lists what is bound and one "Add binding…" row: the prop is picked there,
  // then the … beside it opens the picker on the tray's source
  await page.locator('.u2-designer-properties [data-u2-prop="add-binding"] select').first().selectOption('table');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="add-binding"]').first().click();
  await pickLeaf(page, 'orders (');
  await page.waitForTimeout(600);

  await dropControl(page, 'u2-choice-input', 'choice', 'detailForm');
  await waitStatus(page, 'choiceInput1');
  await selectRow(page, 'choiceInput1');
  await openSection(page, 'Bindings');
  const items = bindField(page, 'items');
  await items.fill('$.orders.columns');
  await items.press('Enter');
  await page.waitForTimeout(400);
  const walked = await bindThroughPicker(page, 'choiceInput1', 'xColumnName', 'plot');
  const bound = findNode((await dump(page)).root, 'choiceInput1')?.bind?.value;
  await shot(page, 'viewers-1-picker');

  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  await page.locator('.u2-designer-surface [data-u2-name="choiceInput1"] select').first().selectOption('total');
  await page.waitForTimeout(800);
  const plot = await viewerAt(page, 'plot', 'xColumnName');
  await shot(page, 'viewers-1-bound-plot');
  await toMode(page, 'Design');
  await page.waitForTimeout(500);
  ok('viewers/1g/a-viewer-look-prop-binds-two-way-through-the-picker',
    walked.labels.includes('plot (u2-viewer-scatter-plot)') && walked.labels.some((l) => l.startsWith('table')) &&
    walked.found === 'xColumnName : string ⇄' && bound === '$.plot.xColumnName' &&
    plot?.type === 'Scatter plot' && plot.props?.xColumnName === 'total' && await specErrors(page) === 0,
    `picked=${walked.found} bind=${bound} plot=${JSON.stringify(plot)} ` +
    `picker=${JSON.stringify(walked.labels.slice(0, 14))}`);
}

/** `leakReport()` (u2/dg, re-exported by U2Demo for this check) plus a census of the live widgets by
 * class and DOM attachment — what names the survivors when the counts disagree. */
const leak = (page) => page.evaluate(() => {
  const pkg = DG.Func.find({package: 'U2demo'})[0]?.package;
  const census = {};
  for (const w of DG.Widget.getAll()) {
    const key = `${w.constructor.name}:${(w.root?.className ?? '').slice(0, 50)}:${w.root?.isConnected ? 'dom' : 'off'}`;
    census[key] = (census[key] ?? 0) + 1;
  }
  return {...(pkg?.getModule?.('package.js') ?? window.u2demo).leakReport(), census};
});

const censusDiff = (a, b) => Object.fromEntries(Object.keys({...a, ...b})
  .filter((k) => (a[k] ?? 0) !== (b[k] ?? 0)).map((k) => [k, `${a[k] ?? 0} -> ${b[k] ?? 0}`]));

const closeViews = async (page) => {
  await page.evaluate(() => {
    for (const view of [...grok.shell.views])
      view.close();
  });
  await page.waitForTimeout(1500);
};

/** The kill-walk of a view close reaches every adopted viewer (VP-7): u2 releases what it owns —
 * the scopes return to what they were before the app opened and no surviving widget is left in the
 * DOM. The platform's own widget count is evidence, not the condition: histogram filter panes leak
 * off-DOM `Legend`/`d4-legend-list` pairs even on a plain `TableView.close()` (D3, core ticket). */
async function checkViewClose(page) {
  await closeViews(page);
  const before = await leak(page);
  await openApp(page, APP.package, APP.func);
  await page.waitForSelector('.u2-designer', {timeout: 60000});
  await page.waitForTimeout(1500);
  await openSpec(page, SAMPLE);
  await waitStatus(page, 'masterDetail');
  await page.waitForTimeout(1500);
  await toMode(page, 'Run');
  await page.waitForTimeout(1000);
  await toMode(page, 'Design');
  await page.waitForTimeout(500);
  const open = await leak(page);
  await closeViews(page);
  const after = await leak(page);
  const survivors = censusDiff(before.census, after.census);
  artifact('viewers-1h-survivors.json', JSON.stringify(survivors, null, 1));
  const counts = (r) => `{scopes: ${r.liveScopes}, widgets: ${r.liveWidgets}}`;
  const inDom = Object.keys(survivors).filter((k) => k.endsWith(':dom'));
  ok('viewers/1h/closing-the-view-releases-every-scope-and-leaves-no-widget-in-the-dom',
    open.liveWidgets > before.liveWidgets && open.liveScopes > before.liveScopes &&
    after.liveScopes === before.liveScopes && inDom.length === 0,
    `before=${counts(before)} open=${counts(open)} after=${counts(after)} ` +
    `DOM-attached survivors=${inDom.length} off-DOM survivors=${JSON.stringify(survivors)}`);
  await reopenApp(page);
}

export const checks = [
  {id: 'viewers/1 design time: the Viewers pane, sample rows in a real grid, no run, canvas selection and the panel',
    run: checkDesignTime},
  {id: 'viewers/2 run mode: the frame as the bus — grid click, row index, a filter click', run: checkRunMode},
  {id: 'viewers/3 Add filter for column… on u2-viewer-filters (WO-V5)', run: checkAddFilter},
  {id: 'viewers/4 panel edit: allowEdit checkbox → viewer and document', run: checkPanelEdit},
  {id: 'viewers/5 a look prop bound through the picker', run: checkBoundLookProp},
  {id: 'viewers/6 view close: scopes back to baseline, no widget left in the DOM', run: checkViewClose},
];
