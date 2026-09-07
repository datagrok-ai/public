/* WO-9 / WO-10 — the editing loop as a user drives it: palette drag, context actions, ribbon,
   keyboard, then the context panel — a property row, a per-child prop, a plain HTML node and a
   binding. The first four run on the round-trip spec, the rest on EDIT_SPEC. */
import {ok, shot} from '../local.mjs';
import {EDIT_SPEC, ROUND_TRIP, balloons, clearBalloons, drag, expandTree, focusTree, labels, named,
  nodeCount, openSpec, paletteItem, panel, panelField, ribbon, row, selectRow, statusText, surfaceNode,
  waitStatus} from '../lib.mjs';

/** The round-trip spec, freshly opened, and the palette unfiltered — the palette count below
 * reads every item. */
export async function fixture(page) {
  await page.locator('.u2-palette input').first().fill('');
  await page.waitForTimeout(150);
  await openSpec(page, ROUND_TRIP);
  await waitStatus(page, 'roundTrip');
  await expandTree(page);
}

/** Set by the palette check — every later check counts against the spec the fixture opened. */
let baseCount = 0;

async function checkPalette(page) {
  baseCount = await nodeCount(page);
  const filter = page.locator('.u2-palette input').first();
  const items = page.locator('.u2-palette-item:visible');
  const all = await items.count();
  await filter.fill('button');
  await page.waitForTimeout(150);
  const filtered = await items.count();
  ok('editing/1a/palette-lists-and-filters', all > 10 && filtered === 1 &&
    await paletteItem(page, 'u2-button').isVisible(),
  `${all} items, ${filtered} after filtering on "button"`);

  await drag(page, paletteItem(page, 'u2-button'), surfaceNode(page, 'rtName'));
  await expandTree(page);
  const after = await named(page);
  const status = await statusText(page);
  await shot(page, 'editing-1-palette-insert');
  ok('editing/1b/drop-inserts-named-node', after.includes('button1'), JSON.stringify(after));
  ok('editing/1c/tree-follows', (await labels(page)).includes('button1'), JSON.stringify(await labels(page)));
  ok('editing/1d/status-counts-and-marks-modified',
    /button1 · \d+ nodes · modified/.test(status) && await nodeCount(page) === baseCount + 1, status);

  // the drop landed between the form's inputs: the canvas shows it exactly there, not at the bottom
  const ys = await page.evaluate(() => Object.fromEntries(['rtName', 'button1', 'rtFlag'].map((n) =>
    [n, document.querySelector(`.u2-designer-surface [data-u2-name="${n}"]`)
      ?.getBoundingClientRect().top ?? -1])));
  ok('editing/1e/canvas-position-matches-the-tree-index',
    ys.rtName >= 0 && ys.rtName < ys.button1 && ys.button1 < ys.rtFlag, JSON.stringify(ys));
}

async function checkActions(page) {
  // the form renders every child — input or not — as a row in place, so canvas order is spec order
  const order = () => page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface [data-u2-name="rtForm"] [data-u2-name]')].map((e) => e.dataset.u2Name));
  const before = await labels(page);
  const box = await surfaceNode(page, 'rtFlag').boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2, {button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu-item').filter({hasText: /^Move Up$/}).first().click();
  await page.waitForTimeout(300);
  const after = await labels(page);
  await ribbon(page, 'Move Up').click();
  await page.waitForTimeout(300);
  const canvas = await order();
  await shot(page, 'editing-2-actions');
  ok('editing/2a/move-up-reorders-document-and-canvas',
    before.indexOf('button1') < before.indexOf('rtFlag') && after.indexOf('rtFlag') < after.indexOf('button1') &&
    JSON.stringify(canvas) === JSON.stringify(['rtFlag', 'rtName', 'button1']),
    `${JSON.stringify(before)} -> ${JSON.stringify(after)} canvas=${JSON.stringify(canvas)}`);

  await selectRow(page, 'rtName');
  await ribbon(page, 'Delete').click();
  await page.waitForTimeout(300);
  const status = await statusText(page);
  ok('editing/2b/delete-removes-and-selects-parent',
    !(await named(page)).includes('rtName') && /› rtForm · /.test(status),
    `status="${status}" named=${JSON.stringify(await named(page))}`);

  await selectRow(page, 'roundTrip');
  ok('editing/2c/root-refuses-delete', await ribbon(page, 'Delete').isDisabled(),
    `Delete disabled=${await ribbon(page, 'Delete').isDisabled()}`);
}

async function checkUndoRedo(page) {
  await focusTree(page);
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  const restored = await named(page);
  ok('editing/3a/ctrl-z-restores', restored.includes('rtName') &&
    (await labels(page)).includes('rtName') && await nodeCount(page) === baseCount + 1,
  `named=${JSON.stringify(restored)} status="${await statusText(page)}"`);

  await ribbon(page, 'Redo').click();
  await page.waitForTimeout(300);
  ok('editing/3b/ribbon-redo-re-deletes', !(await named(page)).includes('rtName'),
    JSON.stringify(await named(page)));

  // the ribbon click parked focus on its button: Ctrl+Z pressed right away — no click into the
  // canvas or the tree — must still be the designer's undo, not the platform's "Nothing to undo"
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  ok('editing/3b2/ctrl-z-lands-straight-after-a-ribbon-click', (await named(page)).includes('rtName'),
    JSON.stringify(await named(page)));
  await ribbon(page, 'Redo').click();
  await page.waitForTimeout(300);

  await focusTree(page);
  for (let i = 0; i < 12 && !await ribbon(page, 'Undo').isDisabled(); i++) {
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(200);
  }
  const status = await statusText(page);
  await shot(page, 'editing-3-undo-redo');
  ok('editing/3c/undo-bottom-is-the-spec-as-opened', await ribbon(page, 'Undo').isDisabled() &&
    !(await named(page)).includes('button1') && await nodeCount(page) === baseCount &&
    !/modified/.test(status), `status="${status}" named=${JSON.stringify(await named(page))}`);
}

/** The affordance is CSS scoped under `.u2-designer:not(.u2-designer-running)`, so run mode drops
 * the styling and keeps the class — the border is what a user sees, and what is asserted. */
async function checkEmptyAffordance(page) {
  await page.locator('.u2-palette input').first().fill('panel');
  await page.waitForTimeout(150);
  await drag(page, paletteItem(page, 'u2-panel'), surfaceNode(page, 'rtFlag'));
  const style = () => page.evaluate(() => {
    const el = document.querySelector('.u2-designer-surface [data-u2-name="panel1"]');
    if (!el)
      return null;
    const css = getComputedStyle(el);
    return {empty: el.classList.contains('u2-designer-empty'), border: css.borderTopStyle,
      height: Math.round(el.getBoundingClientRect().height)};
  });
  const design = await style();
  await shot(page, 'editing-4-empty-affordance');
  ok('editing/4a/empty-container-advertises-itself',
    design !== null && design.empty && design.border === 'dashed' && design.height >= 32,
    JSON.stringify(design));

  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(300);
  const running = await style();
  ok('editing/4b/run-mode-drops-it', running !== null && running.border !== 'dashed', JSON.stringify(running));

  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
  const back = await style();
  ok('editing/4c/back-in-design-mode', back !== null && back.border === 'dashed', JSON.stringify(back));
}

const surfaceText = (page, name) => page.evaluate((name) => document.querySelector(
  `.u2-designer-surface [data-u2-name="${name}"]`)?.textContent.trim() ?? '', name);

async function checkPropertyEdit(page) {
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
  await selectRow(page, 'saveButton');
  const shown = await panel(page);
  ok('editing/5a/panel-offers-editors', shown.Properties?.text === 'Save' && shown.Events?.click === 'cmd:save',
    JSON.stringify(shown));

  const text = panelField(page, 'text');
  await text.click();
  await text.press('End');
  await page.keyboard.type(' now', {delay: 40});
  await page.waitForTimeout(300);
  const live = await page.evaluate(() => ({
    canvas: document.querySelector('.u2-designer-surface [data-u2-name="saveButton"]')?.textContent.trim(),
    inRow: !!document.activeElement?.closest('[data-u2-prop]'),
    typed: document.activeElement?.value,
  }));
  await shot(page, 'editing-5-property-edit');
  ok('editing/5b/edit-reaches-the-canvas-per-keystroke', live.canvas === 'Save now', JSON.stringify(live));
  ok('editing/5c/focus-stays-in-the-row', live.inRow && live.typed === 'Save now', JSON.stringify(live));

  // the row that is already selected: the focus goes to the toolbox without changing what is shown
  await focusTree(page, 'saveButton');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(400);
  const restored = await surfaceText(page, 'saveButton');
  ok('editing/5d/one-undo-restores-the-whole-edit', restored === 'Save', `"${restored}"`);

  // the acceptance repro: after an edit and its undo, the command must still fire in Run mode
  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(300);
  await clearBalloons(page);
  await page.locator('.u2-designer-surface [data-u2-name="saveButton"]').first().click();
  await page.waitForFunction(() => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => /cmd:save ran/.test(b.textContent)), null, {timeout: 5000}).catch(() => {});
  const fired = await balloons(page);
  ok('editing/5e/command-fires-after-edit-and-undo', fired.some((b) => /cmd:save ran/.test(b)),
    JSON.stringify(fired));
  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
}

async function checkChildAndHtml(page) {
  await selectRow(page, 'firstPane');
  const title = panelField(page, 'title');
  await title.fill('Renamed');
  await page.waitForTimeout(400);
  const tabs = await page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface .u2-tabs-label')].map((l) => l.textContent));
  ok('editing/6a/per-child-prop-rebuilds-the-parent', tabs[0] === 'Renamed', JSON.stringify(tabs));

  await selectRow(page, 'heading');
  const shown = await panel(page);
  ok('editing/6b/html-node-has-a-property-model',
    shown.Properties?.text === 'Editing' && 'cls' in (shown.Properties ?? {}),
    JSON.stringify(shown.Properties ?? null));

  await panelField(page, 'text').fill('Edited heading');
  await page.waitForTimeout(400);
  const heading = await surfaceText(page, 'heading');
  await shot(page, 'editing-6-child-and-html');
  ok('editing/6c/html-node-edits-live', heading === 'Edited heading', `"${heading}"`);
}

async function checkBindAndRoundTrip(page) {
  await selectRow(page, 'nameInput');
  // a binding is a guarded field: it commits on Enter/blur, not per keystroke
  await panelField(page, 'value').fill('$.nowhere');
  await panelField(page, 'value').press('Enter');
  await page.waitForTimeout(400);
  const broken = await statusText(page);
  const placeholder = await page.evaluate(() =>
    !!document.querySelector('.u2-designer-surface .u2-spec-error'));
  ok('editing/7a/a-bind-to-nowhere-is-contained', /broken/.test(broken) && placeholder,
    `status="${broken}" placeholder=${placeholder}`);

  await focusTree(page, 'nameInput');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(400);
  const healed = await statusText(page);
  const value = await page.evaluate(() => document.querySelector(
    '.u2-designer-surface [data-u2-name="nameInput"] input')?.value);
  // the bound value is whatever the demo context holds by now — run-mode typing may have appended to it
  ok('editing/7b/undo-brings-the-component-back', !/broken/.test(healed) && /^Aspirin/.test(value ?? ''),
    `status="${healed}" value=${value}`);

  const dump = await openSpec(page);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
  const after = await page.evaluate(() => ({
    named: [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name),
    heading: document.querySelector('.u2-designer-surface [data-u2-name="heading"]')?.textContent,
    tab: document.querySelector('.u2-designer-surface .u2-tabs-label')?.textContent,
    bound: document.querySelector('.u2-designer-surface [data-u2-name="nameInput"] input')?.value,
  }));
  await shot(page, 'editing-7-bind-and-round-trip');
  ok('editing/7c/the-dump-carries-every-edit-and-reopens-identically',
    /"Renamed"/.test(dump) && /Edited heading/.test(dump) && /\$\.reagent/.test(dump) &&
    after.heading === 'Edited heading' && after.tab === 'Renamed' && /^Aspirin/.test(after.bound ?? '') &&
    after.named.includes('saveButton'), `${JSON.stringify(after)} dump=${dump.length} chars`);
}

/** D-8's promised gesture: a double-click on the row label starts the same inline rename F2 does.
 * Runs on the reopened EDIT_SPEC the round-trip check left up. */
async function checkRenameDblclick(page) {
  await selectRow(page, 'secondPane');
  await row(page, 'secondPane').locator('.u2-tree-label').dblclick();
  const editor = page.locator('.u2-tree-rename');
  const open = await editor.count() === 1 && await editor.inputValue() === 'secondPane';
  ok('editing/8a/dblclick-opens-inline-rename', open,
    `editor=${await editor.count()} value="${open ? await editor.inputValue() : ''}"`);

  await editor.fill('sparePane');
  await editor.press('Enter');
  await page.waitForTimeout(400);
  await shot(page, 'editing-8-dblclick-rename');
  // the pane sits in an inactive tab, so it has no canvas element — the tree and the status
  // path are where the rename shows
  ok('editing/8b/rename-commits-through-the-patch', (await labels(page)).includes('sparePane') &&
    /sparePane/.test(await statusText(page)),
  `labels=${JSON.stringify(await labels(page))} status="${await statusText(page)}"`);
}

export const checks = [
  {id: 'editing/1 palette and insert', run: checkPalette},
  {id: 'editing/2 move and delete', run: checkActions},
  {id: 'editing/3 undo and redo', run: checkUndoRedo},
  {id: 'editing/4 empty-container affordance', run: checkEmptyAffordance},
  {id: 'editing/5 property editing in the panel', run: checkPropertyEdit},
  {id: 'editing/6 per-child and HTML props', run: checkChildAndHtml},
  {id: 'editing/7 binding and dump round-trip', run: checkBindAndRoundTrip},
  {id: 'editing/8 rename via double-click', run: checkRenameDblclick},
];
