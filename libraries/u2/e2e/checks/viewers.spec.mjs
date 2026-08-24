/* V-1/V-2/V-2.5 — platform viewers as u2 controls (VIEWERS_AND_WIDGETS.md, viewers/plan.md WO-V4,
   plan-v25.md WO-V12): the master–detail sample with a real DG.Grid, FilterGroup and Form over a
   `sample`-policy source. Design time shows the sample rows without a run; in Run mode the frame
   is the bus — a grid click moves a bound field, a filter click narrows the grid and the form
   together; the panel of a viewer node is the platform property grid over its live look (VP-21),
   where a checkbox click or a column pick is one patch and one undo entry; a look prop binds
   two-way through the picker; a viewer may bind to a control declared after it (R-a); Run-mode
   edits never reach the dump (R-c); `props.filters` rebuilds the panes (VP-31); the fluent factories
   answer the cached wrapper (VP-32); a view close kills every viewer the spec built. */
import {artifact, ok, openApp, shot} from '../local.mjs';
import {APP, balloons, bindField, bindThroughPicker, clearBalloons, dropControl, dumpViaDialog, expandTree,
  fieldValue, focusTree, gridCellCenter, lookField, lookRow, openSection, openSpec, panelField, reopenApp,
  row, selectRow, setField, specErrors, spyRuns, statusText, surfaceNode, toMode, viewerAt, waitCaption,
  waitStatus} from '../lib.mjs';

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

/** The sample with a forward reference (R-a, VP-26): a scatter plot in the detail pane binds its
 * x axis to `xPick`, a choice input over the frame's columns declared after it — in document order
 * the plot comes first, which V-2 refused. */
const FORWARD = (() => {
  const spec = JSON.parse(SAMPLE);
  const detailPane = spec.root.children[2];
  detailPane.children.splice(1, 0, {tag: 'u2-viewer-scatter-plot', name: 'plot',
    bind: {table: '$.orders', xColumnName: '$.xPick'}});
  detailPane.children[2].children.push({tag: 'u2-choice-input', name: 'xPick',
    props: {label: 'X axis', value: 'total'}, bind: {items: '$.orders.columns'}});
  return JSON.stringify(spec, null, 2);
})();

/** The sample with a second frame source: a dropped viewer has two tables to choose from, so the
 * drop seeds no `table` and the panel says what to bind (U1a). */
const TWO_SOURCES = (() => {
  const spec = JSON.parse(SAMPLE);
  spec.components.push({tag: 'u2-func-source', name: 'returns',
    props: {func: 'demoOrders', params: {days: 7}, designData: 'sample', debounce: 50,
      sample: [{orderId: 2001, customer: 'Bayer', city: 'Lviv', total: 90}]}});
  return JSON.stringify(spec, null, 2);
})();

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

const docProps = async (page, name) => findNode((await dump(page)).root, name)?.props;

const docProp = async (page, name, prop) => (await docProps(page, name))?.[prop];

/** The node's whole `props` against what exactly one patch leaves — a document without the key
 * dumps no `props` at all. */
const sameProps = (a, b) => JSON.stringify(a ?? {}) === JSON.stringify(b ?? {});

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

/** The centre of the cell reading `value` in the first d4 grid behind `selector` that has one — the
 * categorical filter pane and the column picker's popup are both canvas, so both are clicked by
 * cell bounds. Hidden columns (no width) are skipped. */
const cellCenterIn = (page, selector, value) => page.evaluate(({selector, value}) => {
  for (const el of document.querySelectorAll(selector)) {
    const grid = DG.Widget.find(el);
    if (!grid?.dataFrame)
      continue;
    for (const column of grid.dataFrame.columns.names()) {
      for (let i = 0; i < grid.dataFrame.rowCount; i++) {
        const cell = grid.cell(column, i);
        if (cell.cell.value !== value)
          continue;
        const b = cell.bounds;
        if (!(b.width > 0))
          break;
        const r = grid.root.getBoundingClientRect();
        return {x: r.left + b.x + b.width / 2, y: r.top + b.y + b.height / 2, column, row: i};
      }
    }
  }
  return null;
}, {selector, value});

const filterCategoryCenter = (page, name, value) =>
  cellCenterIn(page, `.u2-designer-surface [data-u2-name="${name}"] .d4-filter-categorical`, value);

const filterPanes = (page) => page.evaluate(() =>
  document.querySelectorAll('.u2-designer-surface [data-u2-name="filters"] .d4-filter-element').length);

/** What the viewer panel is made of (VP-21): the platform grid's categories and stamped rows under
 * the nearest `Properties` heading before it (the Run-mode hint sits between them, M3), the u2 field a retired form would have stamped, and the wiring panes. */
const lookPanel = (page) => page.evaluate(() => {
  const props = document.querySelector('.u2-designer-properties');
  const grid = props?.querySelector('.u2-designer-look');
  return {
    heading: (() => {
      let el = grid?.previousElementSibling ?? null;
      while (el !== null && el.tagName !== 'H3')
        el = el.previousElementSibling;
      return el?.textContent ?? null;
    })(),
    categories: [...grid?.querySelectorAll('tr.property-grid-category[aria-label]') ?? []]
      .map((r) => r.getAttribute('aria-label')),
    rows: grid?.querySelectorAll('tr.property-grid-item[data-prop-name]').length ?? 0,
    u2Fields: props?.querySelectorAll('[data-u2-prop="allowEdit"]').length ?? 0,
    panes: [...props?.querySelectorAll('.u2-accordion-title') ?? []].map((t) => t.textContent),
  };
});

/** A category a user folded stays folded on the platform grid: the row is revealed the way a user
 * does it — a click on its category row — before its editor is reached for. */
async function revealLookRow(page, name) {
  await page.evaluate((name) => {
    const body = document.querySelector(`.u2-designer-properties .u2-designer-look tr[data-prop-name="${name}"]`)
      ?.closest('tbody');
    if (body?.classList.contains('property-grid-category-body-hide'))
      body.previousElementSibling?.querySelector('tr.property-grid-category')?.click();
  }, name);
  await lookRow(page, name).scrollIntoViewIfNeeded();
}

/** The frame's columns and the look's current pick — what a keyboard step moves between. */
const columnsAt = (page, name, prop) => page.evaluate(({name, prop}) => {
  const viewer = DG.Widget.find(document.querySelector(`.u2-designer-surface [data-u2-name="${name}"]`));
  return {columns: viewer.dataFrame.columns.names(), current: viewer.props[prop] ?? ''};
}, {name, prop});

/** Identity across a re-render: a mark on the live wrapper is gone once the node built a new one. */
const markViewer = (page, name) => page.evaluate((name) => {
  DG.Widget.find(document.querySelector(`.u2-designer-surface [data-u2-name="${name}"]`)).__u2Mark = true;
}, name);

const markedViewer = (page, name) => page.evaluate((name) =>
  DG.Widget.find(document.querySelector(`.u2-designer-surface [data-u2-name="${name}"]`))?.__u2Mark === true, name);

const undo = async (page, name) => {
  await focusTree(page, name);
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(800);
};

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
  const shown = await lookPanel(page);
  const allowEdit = lookField(page, 'allowEdit');
  const unchecked = await allowEdit.count() === 1 && !(await allowEdit.isChecked());
  await shot(page, 'viewers-1-grid-selected');
  ok('viewers/1b/a-canvas-click-selects-the-grid-and-the-panel-is-the-platform-grid-over-its-look',
    /grid · \d+ nodes/.test(status) && unchecked && shown.u2Fields === 0 &&
    ['Data', 'Columns', 'Rows', 'Style'].every((c) => shown.categories.includes(c)) &&
    ['Bindings', 'Events'].every((p) => shown.panes.includes(p)),
    `status="${status}" heading=${shown.heading} categories=${JSON.stringify(shown.categories)} rows=${shown.rows} ` +
    `allowEdit unchecked=${unchecked} u2 fields=${shown.u2Fields} panes=${JSON.stringify(shown.panes)}`);
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

/** VP-23 — a commit in the platform grid is one `set-prop` (one undo entry) that re-creates the
 * viewer; the undo restores the document (the sample holds `allowEdit: false` literally) and re-points
 * the grid onto the restored viewer's look, so its checkbox follows. */
async function checkPanelEdit(page) {
  await selectRow(page, 'grid');
  await revealLookRow(page, 'allowEdit');
  const held = await docProps(page, 'grid');
  await lookField(page, 'allowEdit').click();
  await page.waitForTimeout(800);
  const live = await viewerAt(page, 'grid', 'allowEdit');
  const written = await docProps(page, 'grid');
  const status = await statusText(page);
  await shot(page, 'viewers-1-allow-edit');
  await undo(page, 'grid');
  const undone = await docProps(page, 'grid');
  const restored = await viewerAt(page, 'grid', 'allowEdit');
  const box = await lookField(page, 'allowEdit').isChecked().catch(() => null);
  ok('viewers/1f/a-grid-checkbox-commit-is-one-patch-and-undo-restores-the-look-and-re-points-the-grid',
    live?.props?.allowEdit === true && sameProps(written, {...held, allowEdit: true}) && /modified/.test(status) &&
    sameProps(undone, held) && restored?.props?.allowEdit === false && box === false,
    `document props=${JSON.stringify(held)}; viewer.props.allowEdit=${live?.props?.allowEdit} dump props=` +
    `${JSON.stringify(written)} status="${status}"; after undo: dump props=${JSON.stringify(undone)} ` +
    `viewer=${restored?.props?.allowEdit} checkbox=${box}`);
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
  // the drop bound `table` to the sample's only source (U1b, 1w): the plot is ready for a look bind
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

/** Q2 — the column picker is the grid's one canvas-backed editor: the keyboard steps columns
 * without the popup, the popup is a d4 grid (the `.d4-grid` inside `.d4-column-grid` — the outer element is
 * the `ColumnGrid` wrapper, no frame) clicked by cell bounds; each pick is one patch and
 * one undo entry. The step's direction is read off the live state — `_prev` answers null for the
 * first column (column_combo_box.dart:163-165), so a plot whose auto-pick landed there steps
 * forward instead. */
async function checkColumnPicker(page) {
  await selectRow(page, 'plot');
  await revealLookRow(page, 'xColumnName');
  const held = await docProps(page, 'plot');
  const before = await columnsAt(page, 'plot', 'xColumnName');
  const idx = before.columns.indexOf(before.current);
  const key = idx > 0 ? 'ArrowUp' : 'ArrowDown';
  const expected = before.columns[idx > 0 ? idx - 1 : idx + 1];
  const selector = lookRow(page, 'xColumnName').locator('.d4-column-selector').first();
  await selector.focus();
  await selector.press(key);
  await page.waitForTimeout(800);
  const stepped = await viewerAt(page, 'plot', 'xColumnName');
  const steppedDoc = await docProps(page, 'plot');
  await undo(page, 'plot');
  const undoneStep = await docProps(page, 'plot');
  const afterStep = await viewerAt(page, 'plot', 'xColumnName');
  const keyboard = stepped?.props?.xColumnName === expected && sameProps(steppedDoc, {...held, xColumnName: expected}) &&
    sameProps(undoneStep, held) && afterStep?.props?.xColumnName === before.current;

  await selectRow(page, 'plot');
  await revealLookRow(page, 'xColumnName');
  const now = await columnsAt(page, 'plot', 'xColumnName');
  const target = now.columns.filter((c) => c !== now.current).pop();
  await lookRow(page, 'xColumnName').locator('.d4-column-selector-column').first().click();
  const opened = await page.waitForSelector('.d4-column-grid', {timeout: 3000}).then(() => true, () => false);
  // attached is not painted: the popup grid hit-tests only once it has rendered a frame
  await page.waitForTimeout(500);
  const cell = opened ? await cellCenterIn(page, '.d4-column-grid .d4-grid', target) : null;
  if (cell !== null) {
    await page.mouse.move(cell.x, cell.y);
    await page.waitForTimeout(100);
    await page.mouse.click(cell.x, cell.y);
  }
  await page.waitForTimeout(800);
  const closed = await page.locator('.d4-column-grid').count() === 0;
  const picked = await viewerAt(page, 'plot', 'xColumnName');
  const pickedDoc = await docProps(page, 'plot');
  await shot(page, 'viewers-1-column-pick');
  // a popup left open would eat the Ctrl+Z
  if (!closed)
    await page.keyboard.press('Escape');
  await undo(page, 'plot');
  const undonePick = await docProps(page, 'plot');
  const afterPick = await viewerAt(page, 'plot', 'xColumnName');
  const popup = opened && cell !== null && closed && picked?.props?.xColumnName === target &&
    sameProps(pickedDoc, {...held, xColumnName: target}) && sameProps(undonePick, held) &&
    afterPick?.props?.xColumnName === now.current;
  ok('viewers/1i/a-column-pick-in-the-grid-by-keyboard-and-by-popup-is-one-patch-each',
    keyboard && popup,
    `columns=${JSON.stringify(before.columns)} doc props=${JSON.stringify(held)}; ${key} from "${before.current}": ` +
    `viewer=${stepped?.props?.xColumnName} doc=${JSON.stringify(steppedDoc)} expected=${expected}, undo → doc=` +
    `${JSON.stringify(undoneStep)} viewer=${afterStep?.props?.xColumnName}; popup: opened=${opened} cell=` +
    `${JSON.stringify(cell)} closed=${closed} target=${target} viewer=${picked?.props?.xColumnName} doc=` +
    `${JSON.stringify(pickedDoc)}, undo → doc=${JSON.stringify(undonePick)} viewer=${afterPick?.props?.xColumnName}`);
}

/** R-a (VP-26) — a viewer declared before the control it binds to renders: no placeholder, no
 * "not built" warning, and it follows the control in Run mode; a Design-mode patch on the control
 * re-renders the viewer (a dependent of a named visual node) and it still follows. */
async function checkForwardReference(page) {
  const warnings = [];
  const onConsole = (m) => {
    if (m.type() === 'warning' && /not built/.test(m.text()))
      warnings.push(m.text().slice(0, 200));
  };
  page.on('console', onConsole);
  await openSpec(page, FORWARD);
  await waitStatus(page, 'masterDetail');
  await page.waitForTimeout(1500);
  await expandTree(page);
  page.off('console', onConsole);
  const broken = await specErrors(page);
  const design = await viewerAt(page, 'plot', 'xColumnName');
  await shot(page, 'viewers-1-forward');

  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  await page.locator('.u2-designer-surface [data-u2-name="xPick"] select').first().selectOption('orderId');
  await page.waitForTimeout(800);
  const run = await viewerAt(page, 'plot', 'xColumnName');
  await toMode(page, 'Design');
  await page.waitForTimeout(500);

  await markViewer(page, 'plot');
  await selectRow(page, 'xPick');
  const label = panelField(page, 'label');
  await label.fill('Axis');
  await label.press('Enter');
  await page.waitForTimeout(800);
  const stillMarked = await markedViewer(page, 'plot');
  const after = await viewerAt(page, 'plot', 'xColumnName');
  const written = await docProp(page, 'xPick', 'label');
  ok('viewers/1j/a-viewer-bound-to-a-later-declared-control-renders-and-keeps-following-it',
    broken === 0 && warnings.length === 0 && design?.type === 'Scatter plot' && design.props?.xColumnName === 'total' &&
    run?.props?.xColumnName === 'orderId' && written === 'Axis' && stillMarked === false &&
    after?.props?.xColumnName === 'total',
    `broken=${broken} warnings=${JSON.stringify(warnings)} design=${JSON.stringify(design)} run x=` +
    `${run?.props?.xColumnName}; label → "${written}": re-rendered=${!stillMarked} x=${after?.props?.xColumnName}`);
}

/** R-c (VP-27) — Run-mode interactions change the live controls and nothing else: the dump is the
 * document, so a bound field stays bound and a literal stays what Design mode last wrote. */
async function checkRunModeEdits(page) {
  const held = await docProp(page, 'xPick', 'value');
  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  await setField(page, 'rowIdx', '0');
  await setField(page, 'customerField', 'Zed');
  await page.locator('.u2-designer-surface [data-u2-name="xPick"] select').first().selectOption('orderId');
  await page.waitForTimeout(800);
  const typed = await fieldValue(page, 'customerField');
  const live = await viewerAt(page, 'plot', 'xColumnName');
  await toMode(page, 'Design');
  await page.waitForTimeout(500);
  const doc = (await dump(page)).root;
  const xPick = findNode(doc, 'xPick');
  const customer = findNode(doc, 'customerField');
  const back = await viewerAt(page, 'plot', 'xColumnName');
  ok('viewers/1k/run-mode-edits-never-reach-the-dump',
    held === 'total' && typed === 'Zed' && live?.props?.xColumnName === 'orderId' && xPick?.props?.value === held &&
    customer?.props?.value === undefined && customer?.bind?.value === '$.orders.currentRow.customer' &&
    back?.props?.xColumnName === held,
    `document value=${held}; run: customer="${typed}" x=${live?.props?.xColumnName}; design: dump xPick.value=` +
    `${xPick?.props?.value} customerField.props.value=${customer?.props?.value} bind=${customer?.bind?.value} ` +
    `plot x=${back?.props?.xColumnName}`);
}

const writeFilters = (page, filters) => page.evaluate((filters) => {
  const fg = DG.Widget.find(document.querySelector('.u2-designer-surface [data-u2-name="filters"]'));
  fg.props.filters = filters;
  return fg.getOptions(true).look.filters?.length ?? -1;
}, filters);

/** Q7/VP-31 — `props.filters` written on the live group removes the panes the list no longer
 * names and adds the new ones; the frame's filter is untouched, so the grid and the form keep
 * their counts; no grid is attached to `filters` (the panel shows another node), so the write
 * captures no patch. */
async function checkFiltersWrite(page) {
  const city = {type: 'categorical', column: 'city'};
  const total = {type: 'histogram', column: 'total'};
  const gridBefore = await viewerAt(page, 'grid');
  const panesBefore = await filterPanes(page);
  const status = await statusText(page);
  const oneLook = await writeFilters(page, [city]);
  await page.waitForTimeout(500);
  const one = await filterPanes(page);
  const twoLook = await writeFilters(page, [city, total]);
  await page.waitForTimeout(500);
  const two = await filterPanes(page);
  const grid = await viewerAt(page, 'grid');
  const form = await viewerAt(page, 'form');
  const doc = await docProp(page, 'filters', 'filters');
  await shot(page, 'viewers-1-filters-write');
  ok('viewers/1l/filters-written-through-props-removes-and-adds-panes',
    one === 1 && two === 2 && oneLook === 1 && twoLook === 2 && grid?.filtered === gridBefore?.filtered &&
    form?.filtered === gridBefore?.filtered && doc === undefined && await statusText(page) === status,
    `panes ${panesBefore}→${one}→${two} look.filters ${oneLook}→${twoLook} grid.filtered ${gridBefore?.filtered}→` +
    `${grid?.filtered} form.filtered=${form?.filtered} dump filters=${JSON.stringify(doc)}`);
}

/** VP-32 — `viewers.grid(df)` (exposed by U2Demo as `viewers`, like `leakReport`) answers the
 * cached js-api wrapper, adopted: `DG.Widget.find` and `DG.toJs` meet it, it is a Control, and
 * `viewers.filters` is a real `FilterGroup`. Unhosted viewers are released through the kill-walk,
 * never `close()` (Viewer_Close casts the view). */
const factories = (page) => page.evaluate(() => {
  const pkg = DG.Func.find({package: 'U2demo'})[0]?.package;
  const mod = pkg?.getModule?.('package.js') ?? window.u2demo;
  if (typeof mod?.viewers?.grid !== 'function')
    return {exposed: false};
  const df = grok.data.demo.demog(20);
  const v = mod.viewers.grid(df);
  const fg = mod.viewers.filters(df);
  const out = {exposed: true, found: DG.Widget.find(v.root) === v, toJs: DG.toJs(v.dart) === v,
    control: typeof v.bindStep === 'function', grid: v instanceof DG.Grid, filterGroup: fg instanceof DG.FilterGroup,
    kill: typeof window.grok_Widget_Kill === 'function'};
  for (const w of [v, fg])
    window.grok_Widget_Kill?.(w.root);
  out.gone = DG.Widget.find(v.root) === null && DG.Widget.find(fg.root) === null;
  return out;
});

async function checkFactories(page) {
  const r = await factories(page);
  ok('viewers/1m/the-fluent-factories-adopt-the-cached-wrapper',
    r.exposed && r.found && r.toJs && r.control && r.grid && r.filterGroup && r.kill, JSON.stringify(r));
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
  // a viewer node selected: its panel is a Dart PropertyGrid, killed with the panel (F4)
  await expandTree(page);
  await selectRow(page, 'grid');
  await page.waitForTimeout(500);
  const open = await leak(page);
  await closeViews(page);
  const after = await leak(page);
  const survivors = censusDiff(before.census, after.census);
  artifact('viewers-1h-survivors.json', JSON.stringify(survivors, null, 1));
  const counts = (r) => `{scopes: ${r.liveScopes}, widgets: ${r.liveWidgets}}`;
  const inDom = Object.keys(survivors).filter((k) => k.endsWith(':dom'));
  // evidence, not the condition: the panel's grid is killed on dispose (F4), and what the platform's
  // registry still lists after that is the platform's to answer
  const grids = Object.keys(survivors).filter((k) => k.startsWith('PropertyGrid'));
  ok('viewers/1h/closing-the-view-releases-every-scope-and-leaves-no-widget-in-the-dom',
    open.liveWidgets > before.liveWidgets && open.liveScopes > before.liveScopes &&
    after.liveScopes === before.liveScopes && inDom.length === 0,
    `before=${counts(before)} open=${counts(open)} after=${counts(after)} ` +
    `DOM-attached survivors=${inDom.length} PropertyGrid survivors=${grids.length} off-DOM survivors=${JSON.stringify(survivors)}`);
  await reopenApp(page);
}

/** The UX work order of V-2.5 (progress.md "UX audit (V-2.5)"): the sample opens without filter
 * balloons and with its two panes (B1); the filter group's grid has neither the platform's `table`
 * row nor `columnNames` (M4/M5); "Remove filter…" mirrors Add (M5); a commit keeps focus on its row
 * and Ctrl+Z from the grid's editor undoes it (M2); the rows the document carries are marked and the
 * category rows use the panel's font (taste/P1); Run mode labels the grid live (M3); "Add binding…"
 * reads as words (M8); a viewer dropped without a table opens Bindings with the hint (M6). */
async function checkUxWorkOrder(page) {
  await clearBalloons(page);
  await openSpec(page, SAMPLE);
  await waitStatus(page, 'masterDetail');
  await page.waitForTimeout(1500);
  await expandTree(page);
  const said = await balloons(page);
  const panes = await filterPanes(page);
  ok('viewers/1o/the-sample-opens-without-filter-balloons-and-with-its-two-panes',
    !said.some((t) => /Error adding filter/.test(t)) && panes === 2, `balloons=${JSON.stringify(said)} panes=${panes}`);

  await selectRow(page, 'filters');
  const tableRows = await lookRow(page, 'table').count();
  const aliasRows = await lookRow(page, 'columnNames').count();
  await row(page, 'filters').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  const remove = page.locator('.u2-menu-item').filter({hasText: /^Remove filter…$/});
  const offered = await remove.count() === 1;
  await clearBalloons(page);
  if (offered)
    await remove.first().click();
  else
    await page.keyboard.press('Escape');
  await page.waitForTimeout(500);
  const warned = (await balloons(page)).some((t) => /names no filters/.test(t));
  const noDialog = await page.locator('.u2-dialog').count() === 0;
  ok('viewers/1p/the-filter-group-grid-hides-table-and-columnNames-and-remove-filter-warns-on-an-empty-document',
    tableRows === 0 && aliasRows === 0 && offered && warned && noDialog,
    `table rows=${tableRows} columnNames rows=${aliasRows} offered=${offered} warned=${warned} dialog=${!noDialog}`);

  await row(page, 'filters').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu-item').filter({hasText: /^Add filter for column…$/}).first().click();
  const add = page.locator('.u2-dialog').filter({hasText: 'Column'}).first();
  await add.waitFor({state: 'visible', timeout: 3000});
  await add.locator('select').first().selectOption('orderId');
  await add.locator('button').filter({hasText: /^OK$/i}).first().click();
  await page.waitForTimeout(800);
  const added = (await dump(page)).root.children[0].props?.filters;
  await row(page, 'filters').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu-item').filter({hasText: /^Remove filter…$/}).first().click();
  const dialog = page.locator('.u2-dialog').filter({hasText: 'Filter'}).first();
  const opened = await dialog.waitFor({state: 'visible', timeout: 3000}).then(() => true, () => false);
  let labels = [];
  if (opened) {
    labels = await dialog.locator('select option').allTextContents();
    await dialog.locator('select').first().selectOption({label: 'orderId'});
    await dialog.locator('button').filter({hasText: /^OK$/i}).first().click();
    await page.waitForTimeout(800);
  }
  const removed = (await dump(page)).root.children[0].props?.filters;
  const panesAfter = await filterPanes(page);
  await shot(page, 'viewers-1-remove-filter');
  await undo(page, 'filters');
  await undo(page, 'filters');
  const restored = (await dump(page)).root.children[0].props?.filters;
  ok('viewers/1q/remove-filter-drops-the-picked-entry-from-the-document-and-the-panes-follow',
    opened && added?.length === panes + 1 && labels.includes('orderId') && removed?.length === panes &&
    !removed.some((f) => f.column === 'orderId') && panesAfter === panes && restored === undefined,
    `dialog=${opened} labels=${JSON.stringify(labels)} filters ${JSON.stringify(added)} → ${JSON.stringify(removed)} ` +
    `panes=${panesAfter} after two undos=${JSON.stringify(restored)}`);

  // U6 — the filter group's verbs as buttons under the Node table, for whoever never opens a menu on
  // a tree row: the same action, the same dialog, one patch, one undo entry
  await selectRow(page, 'filters');
  const buttons = page.locator('.u2-designer-properties .u2-designer-verbs button');
  const verbs = await buttons.allTextContents();
  await buttons.filter({hasText: /^Add filter for column…$/}).first().click();
  const addDialog = page.locator('.u2-dialog').filter({hasText: 'Column'}).first();
  const addOpened = await addDialog.waitFor({state: 'visible', timeout: 3000}).then(() => true, () => false);
  if (addOpened) {
    await addDialog.locator('select').first().selectOption('customer');
    await addDialog.locator('button').filter({hasText: /^OK$/i}).first().click();
    await page.waitForTimeout(800);
  }
  const viaButton = (await dump(page)).root.children[0].props?.filters;
  const panesViaButton = await filterPanes(page);
  await shot(page, 'viewers-1-add-filter-button');
  await undo(page, 'filters');
  const panesUndone = await filterPanes(page);
  ok('viewers/1x/the-filter-group-panel-offers-add-and-remove-filter-as-buttons-under-the-node-table',
    verbs.join('|') === 'Add filter for column…|Remove filter…' && addOpened && viaButton?.length === panes + 1 &&
    viaButton[viaButton.length - 1]?.column === 'customer' && panesViaButton === panes + 1 && panesUndone === panes,
    `buttons=${JSON.stringify(verbs)} dialog=${addOpened} filters=${JSON.stringify(viaButton)} panes ${panes}→` +
    `${panesViaButton}→${panesUndone}`);

  await selectRow(page, 'grid');
  await revealLookRow(page, 'allowEdit');
  const styled = await page.evaluate(() => {
    const grid = document.querySelector('.u2-designer-properties .u2-designer-look');
    const category = grid?.querySelector('tr.property-grid-category td');
    const set = grid?.querySelector('tr[data-prop-name="allowEdit"]');
    const other = grid?.querySelector('tr[data-prop-name="showColumnLabels"]');
    return {
      categoryWeight: category ? getComputedStyle(category).fontWeight : null,
      setMarked: set?.classList.contains('u2-designer-look-set') ?? null,
      setWeight: set ? getComputedStyle(set.firstElementChild).fontWeight : null,
      otherMarked: other?.classList.contains('u2-designer-look-set') ?? null,
      nameWidth: set ? set.firstElementChild.getBoundingClientRect().width / set.getBoundingClientRect().width : null,
    };
  });
  ok('viewers/1r/the-rows-the-document-carries-are-bold-and-the-category-rows-use-the-panel-font',
    styled.categoryWeight === '500' && styled.setMarked === true && styled.setWeight === '500' &&
    styled.otherMarked === false && styled.nameWidth !== null && Math.abs(styled.nameWidth - 0.55) < 0.05,
    JSON.stringify(styled));

  const held = await docProps(page, 'grid');
  await lookField(page, 'allowEdit').click();
  await page.waitForTimeout(800);
  const focusedRow = await page.evaluate(() =>
    document.activeElement?.closest('tr[data-prop-name]')?.dataset.propName ?? null);
  const written = await docProps(page, 'grid');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(800);
  const undone = await docProps(page, 'grid');
  const box = await lookField(page, 'allowEdit').isChecked().catch(() => null);
  ok('viewers/1s/a-grid-commit-keeps-focus-on-its-row-and-ctrl-z-from-the-editor-undoes-it',
    focusedRow === 'allowEdit' && sameProps(written, {...held, allowEdit: true}) && sameProps(undone, held) &&
    box === false,
    `focused row after commit=${focusedRow} written=${JSON.stringify(written)} undone=${JSON.stringify(undone)} ` +
    `checkbox=${box}`);

  // U3 — a string commit: the Dart textbox removes itself on Enter before the look is written, which
  // dropped focus to the body and handed the next Ctrl+Z to the platform ("Nothing to undo"); no
  // dialog may open between the Enter and the chord, so the write is read off the live viewer
  await revealLookRow(page, 'title');
  await lookRow(page, 'title').locator('.property-grid-item-view-label').first().click();
  const textbox = lookRow(page, 'title').locator('input').first();
  await textbox.fill('Hello');
  await clearBalloons(page);
  await textbox.press('Enter');
  await page.waitForTimeout(800);
  const typed = await viewerAt(page, 'grid', 'title');
  const focusedAfterEnter = await page.evaluate(() =>
    document.activeElement?.closest('tr[data-prop-name]')?.dataset.propName ?? document.activeElement?.tagName ?? null);
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(800);
  const undoneTitle = await viewerAt(page, 'grid', 'title');
  const complained = await balloons(page);
  const undoneDoc = await docProp(page, 'grid', 'title');
  ok('viewers/1s2/a-string-commit-via-enter-keeps-focus-in-the-grid-and-ctrl-z-from-there-undoes-it',
    typed?.props?.title === 'Hello' && focusedAfterEnter === 'title' && !undoneTitle?.props?.title && undoneDoc === undefined &&
    !complained.some((t) => /Nothing to undo/.test(t)),
    `title written=${typed?.props?.title} focused after Enter=${focusedAfterEnter} after Ctrl+Z: viewer=` +
    `${JSON.stringify(undoneTitle?.props)} doc=${undoneDoc} balloons=${JSON.stringify(complained)}`);

  await openSection(page, 'Bindings');
  const option = await page.locator('.u2-designer-properties [data-u2-prop="add-binding"] option[value="allowEdit"]')
    .first().textContent().catch(() => null);
  ok('viewers/1t/add-binding-lists-the-props-as-words', option === 'Allow Edit', `allowEdit reads "${option}"`);

  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  await selectRow(page, 'grid');
  const live = await lookPanel(page);
  const hint = await page.locator('.u2-designer-properties .u2-designer-hint').filter({hasText: /switch to Design/})
    .count();
  await shot(page, 'viewers-1-run-live-panel');
  await toMode(page, 'Design');
  await page.waitForTimeout(500);
  ok('viewers/1u/run-mode-labels-the-grid-live-with-the-hint', live.heading === 'Properties (live)' && hint === 1,
    `heading=${live.heading} hint=${hint}`);

  // U1b — a viewer dropped into a document with exactly ONE frame source takes its `table` from it
  // in the drop's own patch: it renders bound, and one Ctrl+Z removes the node
  await dropControl(page, 'u2-viewer-bar-chart', 'bar', 'detailPane');
  await waitStatus(page, 'viewerBarChart1');
  await page.waitForTimeout(500);
  const bar = await viewerAt(page, 'viewerBarChart1');
  const seeded = findNode((await dump(page)).root, 'viewerBarChart1')?.bind?.table;
  const broken = await specErrors(page);
  await shot(page, 'viewers-1-bar-auto-bound');
  await undo(page, 'grid');
  const gone = findNode((await dump(page)).root, 'viewerBarChart1') === undefined;
  ok('viewers/1w/a-viewer-dropped-beside-the-only-frame-source-is-bound-to-it-in-the-drop-patch',
    bar?.type === 'Bar chart' && bar.rows === 3 && seeded === '$.orders' && broken === 0 && gone,
    `viewer=${JSON.stringify(bar)} bind.table=${seeded} broken=${broken} one undo removed it=${gone}`);

  // U1a — with two frame sources nothing is seeded; the hint sits at the top of the panel, right
  // under the Node table, never at the end of Bindings
  await openSpec(page, TWO_SOURCES);
  await waitStatus(page, 'masterDetail');
  await page.waitForTimeout(1500);
  await expandTree(page);
  await dropControl(page, 'u2-viewer-bar-chart', 'bar', 'detailPane');
  await waitStatus(page, 'viewerBarChart1');
  await selectRow(page, 'viewerBarChart1');
  const unbound = await page.evaluate(() => {
    const props = document.querySelector('.u2-designer-properties');
    const kids = [...props?.children ?? []];
    const node = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === 'Node');
    const pane = [...props?.querySelectorAll('.u2-accordion-pane') ?? []]
      .find((p) => p.querySelector('.u2-accordion-title')?.textContent === 'Bindings');
    return {
      error: [...props?.querySelectorAll('.u2-table tr') ?? []].some((tr) => tr.children[0]?.textContent === 'Error'),
      expanded: pane?.querySelector('.u2-accordion-header')?.getAttribute('aria-expanded') ?? null,
      hint: kids[node + 2]?.classList.contains('u2-designer-hint') ? kids[node + 2].textContent : null,
      inBindings: pane?.querySelectorAll('.u2-designer-hint').length ?? -1,
    };
  });
  const unseeded = findNode((await dump(page)).root, 'viewerBarChart1')?.bind;
  await shot(page, 'viewers-1-bar-unbound-panel');
  await undo(page, 'grid');
  ok('viewers/1v/a-viewer-dropped-with-two-sources-stays-unbound-and-the-hint-sits-under-the-node-table',
    unbound.error === false && unbound.expanded === 'true' && /Bind `table`/.test(unbound.hint ?? '') &&
    unbound.inBindings === 0 && unseeded === undefined,
    `${JSON.stringify(unbound)} bind=${JSON.stringify(unseeded)}`);
}

export const checks = [
  {id: 'viewers/1 design time: the Viewers pane, sample rows in a real grid, no run, canvas selection and the panel',
    run: checkDesignTime},
  {id: 'viewers/2 run mode: the frame as the bus — grid click, row index, a filter click', run: checkRunMode},
  {id: 'viewers/3 Add filter for column… on u2-viewer-filters (WO-V5)', run: checkAddFilter},
  {id: 'viewers/4 panel edit: allowEdit checkbox in the platform grid → one patch, undo re-points the grid',
    run: checkPanelEdit},
  {id: 'viewers/5 a look prop bound through the picker', run: checkBoundLookProp},
  {id: 'viewers/6 the column picker in the platform grid: keyboard step and popup pick (WO-V12 1i)',
    run: checkColumnPicker},
  {id: 'viewers/7 a forward reference to a later-declared control (R-a, WO-V12 1j)', run: checkForwardReference},
  {id: 'viewers/8 Run-mode edits stay out of the dump (R-c, WO-V12 1k)', run: checkRunModeEdits},
  {id: 'viewers/9 props.filters removes and adds panes (VP-31, WO-V12 1l)', run: checkFiltersWrite},
  {id: 'viewers/10 fluent factories answer the cached wrapper (VP-32, WO-V12 1m)', run: checkFactories},
  {id: 'viewers/12 the UX work order of V-2.5: B1, M2–M6, M8, P1/P2, taste', run: checkUxWorkOrder},
  {id: 'viewers/11 view close: scopes back to baseline, no widget left in the DOM', run: checkViewClose},
];
