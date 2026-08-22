/* P4 WO-17 — the GOAL success criterion, whole: from `New form` to a working master-detail form
   with ZERO hand-written JSON (the wizard, the palette, the binding picker and the function picker
   are the only tools), then that form round-tripped through its own dump, then the design-time
   silence the `schema` policy promises — proven with a function-run spy and a network spy, not by
   looking at the canvas. */
import {artifact, ok, shot} from '../local.mjs';
import {bindThroughPicker, chipMenu, chipNames, dialogOK, dropControl, dumpViaDialog, ensureRecorder,
  expandTree, fieldValue, named, newForm, nodeCount, openSpec, pickFunc, selectChip, selectRow, setField,
  specErrors, spyRuns, toMode, trayMenu, waitStatus} from '../lib.mjs';

export async function fixture(page) {
  await ensureRecorder(page);
  await newForm(page);
}

/** The choice editor of a panel prop — `designData` is a `choices` prop, so it renders a select. */
const panelChoice = (page, name) =>
  page.locator(`.u2-designer-properties [data-u2-prop="${name}"] select`).first();

async function checkFullBuild(page) {
  const blank = await nodeCount(page);
  ok('full-build/1a/new-form-starts-from-blank', blank === 1 &&
    (await chipNames(page)).length === 0 && await specErrors(page) === 0,
  `${blank} node(s), chips=${JSON.stringify(await chipNames(page))}`);

  // 1. the data: the tray wizard, a picked function and a literal parameter
  await trayMenu(page, 'Function or query…');
  await pickFunc(page, 'demoOrders');
  await page.locator('.u2-func-picker [data-u2-prop="days"] input').first().fill('30');
  await dialogOK(page, '.u2-func-picker');
  await page.waitForTimeout(600);
  await shot(page, 'full-build-1-source');
  ok('full-build/1b/the-wizard-puts-the-function-on-the-tray',
    (await chipNames(page)).includes('funcSource1'), JSON.stringify(await chipNames(page)));

  // live + one explicit Refresh: what makes the frame — and so its columns — visible to the picker
  await selectChip(page, 'funcSource1');
  await panelChoice(page, 'designData').selectOption('live');
  await page.waitForTimeout(500);
  await chipMenu(page, 'funcSource1', 'Refresh');

  // 2. the form: four controls, dropped
  await dropControl(page, 'u2-number-input', 'number', 'form1');
  await dropControl(page, 'u2-text-input', 'text', 'numberInput1');
  await dropControl(page, 'u2-text-input', 'text', 'textInput1');
  await dropControl(page, 'u2-button', 'button', 'textInput2');
  await expandTree(page);
  const built = await named(page);
  await shot(page, 'full-build-1-controls');
  ok('full-build/1c/the-palette-builds-the-form', JSON.stringify(built) ===
    JSON.stringify(['form1', 'numberInput1', 'textInput1', 'textInput2', 'button1']),
  JSON.stringify(built));

  // 3. the wiring: three binds through the picker, one event through the function picker
  await bindThroughPicker(page, 'numberInput1', 'currentRowIdx', 'funcSource1');
  await bindThroughPicker(page, 'textInput1', 'customer', 'funcSource1', 'currentRow :');
  const rows = await bindThroughPicker(page, 'textInput2', 'rowCount', 'funcSource1');
  const wired = JSON.parse(await dumpViaDialog(page)).root.children;
  await shot(page, 'full-build-1-bound');
  ok('full-build/1d/the-picker-binds-the-navigator-the-detail-field-and-the-row-count',
    wired[0].bind?.value === '$.funcSource1.currentRowIdx' &&
    wired[1].bind?.value === '$.funcSource1.currentRow.customer' &&
    wired[2].bind?.value === '$.funcSource1.rowCount',
    `${JSON.stringify(wired.map((c) => c.bind))} last picked=${rows.found}`);

  await selectRow(page, 'button1');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="click"]').first().click();
  await pickFunc(page, 'u2Record');
  await page.locator('.u2-func-picker [data-u2-prop="text"] input').first().fill('saved');
  await dialogOK(page, '.u2-func-picker');
  await page.waitForTimeout(500);

  // 4. Run: a master-detail form, and nothing in it was typed as JSON
  await toMode(page, 'Run');
  await page.waitForTimeout(900);
  await setField(page, 'numberInput1', '1');
  const detail = await fieldValue(page, 'textInput1');
  const counted = await fieldValue(page, 'textInput2');
  await page.evaluate(() => window.__u2Fired = '');
  await page.locator('.u2-designer-surface [data-u2-name="button1"]').first().click();
  await page.waitForTimeout(600);
  const fired = await page.evaluate(() => window.__u2Fired);
  await shot(page, 'full-build-1-run');
  ok('full-build/1e/run-mode-is-a-working-master-detail-form-with-no-hand-written-json',
    detail === 'Bayer' && counted === '4' && fired === 'saved' && await specErrors(page) === 0,
    `detail="${detail}" rows="${counted}" fired="${fired}"`);

  // 5. the round trip: the form's own dump reopens to the same form
  await toMode(page, 'Design');
  const dump = await dumpViaDialog(page);
  artifact('full-build-1-dump.json', dump);
  await openSpec(page);
  await waitStatus(page, 'form1');
  await expandTree(page);
  const again = await dumpViaDialog(page);
  await shot(page, 'full-build-1-round-trip');
  ok('full-build/1f/copy-spec-then-open-restores-the-identical-canvas',
    again === dump && JSON.stringify(await named(page)) === JSON.stringify(built) &&
    await specErrors(page) === 0, `${dump.length} chars, reopened ${again.length}`);

  // 6. design-time silence: a policy flip must empty the preview without running anything
  await spyRuns(page);
  const requests = [];
  const record = (request) => requests.push(request.url());
  page.on('request', record);
  await selectChip(page, 'funcSource1');
  await chipMenu(page, 'funcSource1', 'Refresh');
  // the row count, not a field of the current row: a frame the designer has only just built has
  // no current row until something sets one, so `currentRow.*` is honestly empty either way
  const filled = await fieldValue(page, 'textInput2');
  const onRefresh = await page.evaluate(() =>
    window.__u2Runs.filter((name) => name === 'demoOrders').length);

  await page.evaluate(() => window.__u2Runs = []);
  requests.length = 0;
  await panelChoice(page, 'designData').selectOption('schema');
  await page.waitForTimeout(1500);
  // no frame at all under `schema`: the row count falls back to 0 and every field of the current
  // row reads empty — the preview holds nothing a run put there
  const emptied = await fieldValue(page, 'textInput2');
  const detailGone = await fieldValue(page, 'textInput1');
  const onFlip = await page.evaluate(() => window.__u2Runs);
  page.off('request', record);
  await page.evaluate(() => window.__u2Sub?.unsubscribe?.());
  await shot(page, 'full-build-1-schema');
  ok('full-build/1g/the-schema-policy-empties-the-preview-and-runs-nothing',
    filled === '4' && onRefresh >= 1 && emptied === '0' && detailGone === '' &&
    onFlip.length === 0 && requests.length === 0,
    `refreshed="${filled}" (${onRefresh} run(s)) → schema rows="${emptied}" detail="${detailGone}" ` +
    `runs=${JSON.stringify(onFlip)} requests=${JSON.stringify(requests.slice(0, 5))}`);
}

export const checks = [
  {id: 'full-build/1 the whole walk: New form → wizard → palette → pickers → Run → round trip → silence',
    run: checkFullBuild},
];
