/* P4 WO-12/14/16 and the acceptance feedback — the components tray and the data sources: a drag
   that lands on the tray, the chip as a non-visual component's handle, a `u2-state` bound into two
   inputs, a func-source driving a master-detail form, an entity-source over the users fixture, the
   guided ways onto the tray (workspace picker, source wizard, function picker on an event), and
   the source's own feedback (unrun columns, Refresh reporting, out-of-range rows). */
import {ok, shot} from '../local.mjs';
import {balloons, bindField, bindThroughPicker, chipMenu, chipNames, clearBalloons, dialogCancel, dialogOK,
  drag, dropControl, dumpViaDialog, ensureDemog, ensureRecorder, expandTree, fieldValue, focusTree, labels,
  newForm, nodeCount, openSection, openSpec, paletteItem, panel, pickFunc, pickLeaf, pickerLabels, pickerRow,
  pickerWalk, selectChip, selectRow, setField, setPanelField, specErrors, statusText, toMode, trayMenu,
  waitStatus} from '../lib.mjs';

const TRAY_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  root: {tag: 'u2-form', name: 'trayForm', children: [
    {tag: 'u2-text-input', name: 'firstInput', props: {label: 'First'}, bind: {value: '$.nowhere'}},
    {tag: 'u2-text-input', name: 'secondInput', props: {label: 'Second'}, bind: {value: '$.nowhere'}},
  ]},
}, null, 2);

export async function fixture(page) {
  await page.locator('.u2-palette input').first().fill('');
  await page.waitForTimeout(150);
  await openSpec(page, TRAY_SPEC);
  await waitStatus(page, 'trayForm');
  await expandTree(page);
}

async function checkComponentsTray(page) {
  const panes = await page.evaluate(() => [...document.querySelectorAll('.u2-palette .u2-accordion-title')]
    .map((t) => t.textContent));
  ok('tray-and-sources/1a/palette-groups-the-sources-under-data', panes.includes('Data'), JSON.stringify(panes));

  // an unfiltered u2-state sits below the fold
  await page.locator('.u2-palette input').first().fill('state');
  await page.waitForTimeout(200);
  const before = await nodeCount(page);
  await drag(page, paletteItem(page, 'u2-state'), page.locator('.u2-designer-tray').first());
  const status = await waitStatus(page, 'state1');
  await expandTree(page);
  await shot(page, 'tray-and-sources-1-tray-drop');
  ok('tray-and-sources/1b/a-data-drag-lands-on-the-tray', (await chipNames(page)).includes('state1') &&
    await nodeCount(page) === before + 1, `chips=${JSON.stringify(await chipNames(page))} status="${status}"`);
  ok('tray-and-sources/1c/the-tree-gains-the-component-row', (await labels(page)).includes('state1'),
    JSON.stringify(await labels(page)));

  const caption = await selectChip(page, 'state1');
  const shown = await panel(page);
  ok('tray-and-sources/1d/the-chip-selects-the-component-and-the-panel-describes-it',
    caption === 'state1 (u2-state)' && shown.Properties?.type === 'string' &&
    'initial' in (shown.Properties ?? {}),
  `caption="${caption}" props=${JSON.stringify(shown.Properties ?? null)}`);

  const broken = await specErrors(page);
  await selectRow(page, 'firstInput');
  await setPanelField(page, 'value', '$.state1');
  await selectRow(page, 'secondInput');
  await setPanelField(page, 'value', '$.state1');
  await shot(page, 'tray-and-sources-1-bound');
  ok('tray-and-sources/1e/binding-to-a-component-through-the-panel-resolves',
    broken === 2 && await specErrors(page) === 0,
    `broken before=${broken} after=${await specErrors(page)}`);

  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(400);
  const first = page.locator('.u2-designer-surface [data-u2-name="firstInput"] input').first();
  const second = page.locator('.u2-designer-surface [data-u2-name="secondInput"] input').first();
  await first.click();
  await page.keyboard.type('shared', {delay: 20});
  await page.waitForTimeout(200);
  const echoed = await second.inputValue();
  ok('tray-and-sources/1f/the-state-carries-the-value-between-both-inputs-in-run-mode', echoed === 'shared',
    `first="${await first.inputValue()}" second="${echoed}"`);

  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
  await focusTree(page, 'trayForm');
  for (let i = 0; i < 3; i++) {
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(250);
  }
  const undone = await statusText(page);
  await shot(page, 'tray-and-sources-1-undone');
  ok('tray-and-sources/1g/ctrl-z-unwinds-the-binds-and-the-component', (await chipNames(page)).length === 0 &&
    await nodeCount(page) === before && await specErrors(page) === 2,
  `chips=${JSON.stringify(await chipNames(page))} status="${undone}"`);
}

/* The GOAL's master-detail criterion minus the list (no df-consuming list component exists yet,
   Risk R6): a func-source over a demo table, a row navigator bound two-way to its `currentRowIdx`,
   and fields bound to the current row's columns, one of which proves the write reaches the frame
   by echoing into a second field. */
const ORDERS_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-func-source', name: 'orders',
    props: {func: 'demoOrders', params: {days: 30}, designData: 'live', debounce: 50}}],
  root: {tag: 'u2-form', name: 'ordersForm', children: [
    {tag: 'u2-number-input', name: 'rowIdx', props: {label: 'Row', mode: 'int'},
      bind: {value: '$.orders.currentRowIdx'}},
    {tag: 'u2-text-input', name: 'customerField', props: {label: 'Customer'},
      bind: {value: '$.orders.currentRow.customer'}},
    {tag: 'u2-text-input', name: 'customerEcho', props: {label: 'Echo'},
      bind: {value: '$.orders.currentRow.customer'}},
    {tag: 'u2-text-input', name: 'cityField', props: {label: 'City'},
      bind: {value: '$.orders.currentRow.city'}},
  ]},
}, null, 2);

async function checkFuncSource(page) {
  await openSpec(page, ORDERS_SPEC);
  await waitStatus(page, 'ordersForm');
  await page.waitForTimeout(600);
  await shot(page, 'tray-and-sources-2-design');
  ok('tray-and-sources/2a/design-time-silence-for-a-non-query', await fieldValue(page, 'customerField') === '' &&
    await specErrors(page) === 0,
  `customer="${await fieldValue(page, 'customerField')}" broken=${await specErrors(page)}`);

  await toMode(page, 'Run');
  await setField(page, 'rowIdx', '0');
  await shot(page, 'tray-and-sources-2-run');
  ok('tray-and-sources/2b/run-mode-runs-the-source-and-fills-the-bound-fields',
    await fieldValue(page, 'customerField') === 'Aspirin Labs' &&
    await fieldValue(page, 'cityField') === 'Kyiv',
    `customer="${await fieldValue(page, 'customerField')}" city="${await fieldValue(page, 'cityField')}"`);

  await setField(page, 'rowIdx', '1');
  ok('tray-and-sources/2c/stepping-the-row-index-moves-every-bound-field',
    await fieldValue(page, 'customerField') === 'Bayer' &&
    await fieldValue(page, 'cityField') === 'Lviv',
    `customer="${await fieldValue(page, 'customerField')}" city="${await fieldValue(page, 'cityField')}"`);

  await setField(page, 'customerField', 'Bayer AG');
  await shot(page, 'tray-and-sources-2-edited');
  ok('tray-and-sources/2d/editing-a-bound-field-writes-the-frame',
    await fieldValue(page, 'customerEcho') === 'Bayer AG',
    `echo="${await fieldValue(page, 'customerEcho')}"`);

  await setField(page, 'rowIdx', '0');
  await setField(page, 'rowIdx', '1');
  ok('tray-and-sources/2e/the-edit-lives-in-the-frame-not-in-the-field',
    await fieldValue(page, 'customerField') === 'Bayer AG',
    `customer="${await fieldValue(page, 'customerField')}"`);
  await toMode(page, 'Design');
}

/* `u2-entity-source` against the local-mode `users` fixture (web/local/api.json's canned
   responses): a real dapi listing, paged by the platform's own source, feeding a choice input's
   items. Always live, so it works in Design mode too. */
const USERS_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-entity-source', name: 'people',
    props: {entity: 'users', pageSize: 10}}],
  root: {tag: 'u2-form', name: 'usersForm', children: [
    {tag: 'u2-choice-input', name: 'userPick', props: {label: 'User'},
      bind: {items: '$.people.names'}},
    {tag: 'u2-text-input', name: 'stateField', props: {label: 'State'},
      bind: {value: '$.people.state'}},
  ]},
}, null, 2);

async function checkEntitySource(page) {
  await openSpec(page, USERS_SPEC);
  await waitStatus(page, 'usersForm');
  await page.waitForTimeout(1200);
  await shot(page, 'tray-and-sources-3-users');

  const options = await page.evaluate(() =>
    [...document.querySelectorAll('.u2-designer-surface [data-u2-name="userPick"] option')]
      .map((o) => o.textContent));
  ok('tray-and-sources/3a/the-collection-reaches-a-settled-state', await fieldValue(page, 'stateField') === 'done',
    `state="${await fieldValue(page, 'stateField')}" broken=${await specErrors(page)}`);
  ok('tray-and-sources/3b/names-feed-the-choice-inputs-items',
    options.includes('Alice Mendeleev') && options.includes('Bob Curie'),
    JSON.stringify(options));
  ok('tray-and-sources/3c/nothing-is-broken', await specErrors(page) === 0, `broken=${await specErrors(page)}`);
}

/* The guided ways onto the tray and the function picker: the `+` menu's workspace picker, the
   source wizard with a param bound to an input, and an event wired to a platform function with a
   literal argument — none of it hand-written JSON. */
const WIZARD_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  root: {tag: 'u2-form', name: 'wizardForm', children: [
    {tag: 'u2-number-input', name: 'daysInput', props: {label: 'Days', mode: 'int', value: 30}},
    {tag: 'u2-text-input', name: 'siteField', props: {label: 'Site'}},
    {tag: 'u2-text-input', name: 'countField', props: {label: 'Orders'}},
    {tag: 'u2-button', name: 'fireButton', props: {text: 'Fire'}},
  ]},
}, null, 2);

async function checkSourceWizard(page) {
  // what the wizard offers: a table the user already has open, and a function to fire from a button
  await ensureDemog(page);
  await ensureRecorder(page);
  await openSpec(page, WIZARD_SPEC);
  await waitStatus(page, 'wizardForm');
  await expandTree(page);

  const before = await nodeCount(page);
  await trayMenu(page, 'Open tables', 'u2demog');
  await shot(page, 'tray-and-sources-4-table-source');
  ok('tray-and-sources/4a/the-workspace-picker-inserts-a-table-source',
    (await chipNames(page)).includes('u2demog1') && await nodeCount(page) === before + 1,
    `chips=${JSON.stringify(await chipNames(page))}`);

  await selectRow(page, 'siteField');
  await openSection(page, 'Bindings');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="value"]').first().click();
  await page.waitForSelector('.u2-bind-picker', {timeout: 5000});
  await pickerRow(page, 'u2demog1').locator('.u2-tree-twistie').click();
  await page.waitForTimeout(300);
  // 'currentRow :' — `currentRowIdx` is a leaf above it, and it carries the same prefix
  await pickerRow(page, 'currentRow :').locator('.u2-tree-twistie').click();
  await page.waitForTimeout(300);
  // the columns of the open table follow their walkable step, in the frame's own order
  const labels = await pickerLabels(page);
  const column = labels[labels.findIndex((label) => label.startsWith('currentRow :')) + 1] ?? '';
  if (column !== '')
    await pickerRow(page, column).click();
  await dialogOK(page, '.u2-bind-picker');
  await page.waitForTimeout(400);
  const bound = JSON.parse(await dumpViaDialog(page)).root.children[1].bind?.value;
  const shown = await fieldValue(page, 'siteField');
  await shot(page, 'tray-and-sources-4-table-bound');
  ok('tray-and-sources/4b/the-picker-walks-the-open-table-to-a-column-of-its-current-row',
    /^\$\.u2demog1\.currentRow[.[]/.test(bound ?? '') && shown !== '',
    `bind=${bound} value="${shown}" column="${column}" rows=${JSON.stringify(labels.slice(0, 12))}`);

  await trayMenu(page, 'Function or query…');
  await pickFunc(page, 'demoOrders');
  await page.locator('.u2-func-picker [data-u2-bind-pick="days"]').first().click();
  const param = await pickLeaf(page, 'value', 'daysInput');
  await dialogOK(page, '.u2-func-picker');
  await page.waitForTimeout(800);
  const wizard = JSON.parse(await dumpViaDialog(page));
  const source = (wizard.components ?? []).find((c) => c.name === 'funcSource1');
  await shot(page, 'tray-and-sources-4-func-source');
  ok('tray-and-sources/4c/the-wizard-inserts-the-function-with-its-param-bound',
    (await chipNames(page)).includes('funcSource1') && source?.props.func === 'demoOrders' &&
    source?.bind?.['params.days'] === '$.daysInput.value' && !('designData' in (source?.props ?? {})),
    `${JSON.stringify(source)} picked=${param.found}`);

  await selectRow(page, 'countField');
  await openSection(page, 'Bindings');
  const path = bindField(page, 'value');
  await path.fill('$.funcSource1.rowCount');
  await path.press('Enter');
  await page.waitForTimeout(400);
  await toMode(page, 'Run');
  await page.waitForTimeout(1000);
  const wide = await fieldValue(page, 'countField');
  await setField(page, 'daysInput', '5');
  await page.waitForTimeout(1000);
  const narrow = await fieldValue(page, 'countField');
  await shot(page, 'tray-and-sources-4-rerun');
  ok('tray-and-sources/4d/editing-the-bound-input-re-runs-the-source', wide === '4' && narrow === '2',
    `30 days="${wide}" 5 days="${narrow}"`);

  await toMode(page, 'Design');
  await selectRow(page, 'fireButton');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="click"]').first().click();
  await pickFunc(page, 'u2Record');
  await page.locator('.u2-func-picker [data-u2-prop="text"] input').first().fill('fired');
  await dialogOK(page, '.u2-func-picker');
  await page.waitForTimeout(500);
  const wired = JSON.parse(await dumpViaDialog(page)).root.children[3].on?.click;
  await toMode(page, 'Run');
  await page.locator('.u2-designer-surface [data-u2-name="fireButton"]').first().click();
  await page.waitForTimeout(600);
  const fired = await page.evaluate(() => window.__u2Fired);
  await shot(page, 'tray-and-sources-4-event-fired');
  ok('tray-and-sources/4e/the-picker-wires-an-event-to-a-function-and-run-mode-fires-it',
    wired?.cmd === 'cmd:u2Record' && wired?.args?.text === 'fired' && fired === 'fired',
    `${JSON.stringify(wired)} fired="${fired}"`);
  await toMode(page, 'Design');
}

/* P4 acceptance — the source's own feedback. What the pass found: `currentRow` opened onto nothing
   at all ("I genuinely thought my data source was broken"), Refresh changed nothing anywhere, and a
   row index past the last row left the previous row's text sitting in the fields as if it were
   data. Every one of those is a question the designer must answer out loud. */
async function checkSourceFeedback(page) {
  await newForm(page);
  await trayMenu(page, 'Function or query…');
  await pickFunc(page, 'demoOrders');
  await page.locator('.u2-func-picker [data-u2-prop="days"] input').first().fill('30');
  await dialogOK(page, '.u2-func-picker');
  await page.waitForTimeout(600);
  await dropControl(page, 'u2-number-input', 'number', 'form1');
  await dropControl(page, 'u2-text-input', 'text', 'numberInput1');
  await expandTree(page);

  // a package function never auto-runs at design time, so its frame — and its columns — are honestly
  // unknown until something runs it. The branch must SAY that instead of opening onto silence
  await selectRow(page, 'textInput1');
  await openSection(page, 'Bindings');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="value"]').first().click();
  const unrun = await pickerWalk(page, 'funcSource1', 'orders', 'currentRow :');
  await dialogCancel(page, '.u2-bind-picker');
  await page.waitForTimeout(300);
  await shot(page, 'tray-and-sources-5-unrun');
  ok('tray-and-sources/5a/an-unrun-sources-columns-explain-themselves-instead-of-showing-nothing',
    unrun.some((label) => label.startsWith('nothing here yet')) &&
    !unrun.some((label) => label.startsWith('customer')), JSON.stringify(unrun));

  await selectChip(page, 'funcSource1');
  const idle = await panel(page);
  ok('tray-and-sources/5b/the-panel-says-the-source-has-not-run',
    (idle.Status?.State ?? '').startsWith('idle'), JSON.stringify(idle.Status ?? null));

  await clearBalloons(page);
  await chipMenu(page, 'funcSource1', 'Refresh');
  await page.waitForFunction(() => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => /funcSource1/.test(b.textContent)), null, {timeout: 8000}).catch(() => {});
  const said = await balloons(page);
  const ready = await panel(page);
  await shot(page, 'tray-and-sources-5-refreshed');
  ok('tray-and-sources/5c/refresh-runs-the-source-and-names-what-came-back',
    said.some((text) => /funcSource1: ready/.test(text)) && ready.Status?.Rows === '4',
    `balloons=${JSON.stringify(said)} status=${JSON.stringify(ready.Status ?? null)}`);

  const walked = await bindThroughPicker(page, 'textInput1', 'customer',
    'funcSource1', 'orders', 'currentRow :');
  ok('tray-and-sources/5d/a-source-that-has-run-lists-its-columns',
    walked.found !== undefined && !walked.labels.some((l) => l.startsWith('nothing here yet')),
    `found=${walked.found} labels=${JSON.stringify(walked.labels)}`);
  await bindThroughPicker(page, 'numberInput1', 'currentRowIdx', 'funcSource1', 'orders');

  await toMode(page, 'Run');
  await page.waitForTimeout(900);
  await setField(page, 'numberInput1', '1');
  const inRange = await fieldValue(page, 'textInput1');
  await setField(page, 'numberInput1', '9');
  const outside = await fieldValue(page, 'textInput1');
  await setField(page, 'numberInput1', '1');
  const back = await fieldValue(page, 'textInput1');
  await shot(page, 'tray-and-sources-5-out-of-range');
  ok('tray-and-sources/5e/a-row-index-past-the-last-row-reads-empty-both-ways',
    inRange === 'Bayer' && outside === '' && back === 'Bayer',
    `row 1="${inRange}" row 9="${outside}" back="${back}"`);
  await toMode(page, 'Design');
}

export const checks = [
  {id: 'tray-and-sources/1 components tray: Data palette, drag to tray, chip selection, state binding',
    run: checkComponentsTray},
  {id: 'tray-and-sources/2 func-source master-detail: row navigator, bound fields, write-back', run: checkFuncSource},
  {id: 'tray-and-sources/3 entity-source over the users fixture: state, names into a choice input',
    run: checkEntitySource},
  {id: 'tray-and-sources/4 tray + menu: workspace picker, source wizard, function picker on an event',
    run: checkSourceWizard},
  {id: 'tray-and-sources/5 source feedback: the unrun columns note, Refresh reporting, the panel status, ' +
    'out-of-range rows', run: checkSourceFeedback},
];
