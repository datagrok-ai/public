/* P4 WO-15 and the acceptance fixes — the binding picker: the Bindings section as the place a
   binding is ADDED, the … button walking the tray source to its value, the cycle refusal a
   hand-typed path still gets; then both pickers under a query — the match revealed, the expansion
   kept, the empty states, the function row's layout. */
import {ok, shot} from '../local.mjs';
import {balloons, bindField, bindThroughPicker, clearBalloons, dialogCancel, dumpViaDialog, ensureDemog,
  expandTree, focusTree, openSection, openSpec, panel, pickerLabels, pickerRow, selectRow, specErrors,
  waitStatus} from '../lib.mjs';

const PICK_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-state', name: 'draft', props: {type: 'string', initial: 'seeded'}}],
  root: {tag: 'u2-form', name: 'pickForm', children: [
    {tag: 'u2-text-input', name: 'pickInput', props: {label: 'Pick'}},
  ]},
}, null, 2);

export async function fixture(page) {
  await openSpec(page, PICK_SPEC);
  await waitStatus(page, 'pickForm');
  await expandTree(page);
}

const boundPath = (dump, name) => JSON.parse(dump).root.children[0].bind?.[name];

async function checkBindPicker(page) {
  await selectRow(page, 'pickInput');
  const shown = await panel(page);
  ok('pickers/1a/bindings-lists-bound-rows-only-plus-add-binding', shown.Bindings?.value === undefined &&
    'add-binding' in (shown.Bindings ?? {}), JSON.stringify(shown.Bindings ?? null));

  await openSection(page, 'Bindings');
  await page.locator('[data-u2-bind-pick="value"]').first().click();
  await page.waitForSelector('.u2-bind-picker', {timeout: 5000});
  await pickerRow(page, 'draft').locator('.u2-tree-twistie').click();
  await page.waitForTimeout(300);
  const rows = await pickerLabels(page);
  await shot(page, 'pickers-1-picker');
  ok('pickers/1b/the-picker-walks-the-tray-source-to-its-value',
    rows.some((r) => r.startsWith('draft')) && rows.some((r) => r.startsWith('value')),
    JSON.stringify(rows));

  await pickerRow(page, 'value').click();
  await page.locator('.u2-bind-picker button').filter({hasText: /^OK$/}).first().click();
  await page.waitForTimeout(400);
  const bound = await dumpViaDialog(page);
  const value = await page.evaluate(() =>
    document.querySelector('.u2-designer-surface [data-u2-name="pickInput"] input')?.value);
  await shot(page, 'pickers-1-bound');
  ok('pickers/1c/picking-a-leaf-binds-and-the-canvas-follows',
    boundPath(bound, 'value') === '$.draft.value' && value === 'seeded' &&
    await specErrors(page) === 0, `bind=${boundPath(bound, 'value')} value="${value}"`);

  await focusTree(page, 'pickInput');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(400);
  const undone = await dumpViaDialog(page);
  ok('pickers/1d/one-undo-restores-the-unbound-input', boundPath(undone, 'value') === undefined,
    undone.slice(0, 200));

  await clearBalloons(page);
  await bindThroughPicker(page, 'pickInput', 'value', 'draft');
  const field = bindField(page, 'value');
  await field.fill('$.pickInput');
  await field.press('Enter');
  await page.waitForTimeout(400);
  const refused = (await balloons(page)).find((b) => /cycle/.test(b)) ?? '';
  await shot(page, 'pickers-1-cycle-refused');
  ok('pickers/1e/a-self-referential-path-is-refused-with-the-loop',
    refused.includes('pickInput → pickInput') && await field.inputValue() === '$.pickInput' &&
    boundPath(await dumpViaDialog(page), 'value') === '$.draft.value', `balloon="${refused}"`);
  await clearBalloons(page);
}

/* The binding tree used to filter its matches into branches it then left closed, and to reset
   every expansion the moment the box was cleared; the function list drew the platform's entity
   markup, whose icon sat on top of the caption and whose descenders fell outside the row, and
   matched on text no row showed. */
const SEARCH_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-table-source', name: 'demog', props: {table: 'u2demog'}}],
  root: {tag: 'u2-form', name: 'searchForm', children: [
    {tag: 'u2-text-input', name: 'searchField', props: {label: 'Site'}},
    {tag: 'u2-button', name: 'searchButton', props: {text: 'Fire'}},
  ]},
}, null, 2);

const pickerSearch = (page, dialog) => page.locator(`${dialog} input`).first();

/** The geometry the acceptance shot showed going wrong: the icon over the caption, and the caption's
 * descenders below the row's own bottom edge. */
const funcRowLayout = (page) => page.evaluate(() => {
  const row = document.querySelector('.u2-func-picker .u2-list-row');
  const label = row?.querySelector('.u2-func-row-label');
  if (!row || !label)
    return null;
  const r = row.getBoundingClientRect();
  const l = label.getBoundingClientRect();
  const i = row.querySelector('.u2-func-row-icon')?.getBoundingClientRect();
  return {text: label.textContent, gap: i ? Math.round(l.left - i.right) : 0,
    top: Math.round(l.top - r.top), bottom: Math.round(r.bottom - l.bottom),
    description: row.querySelector('.u2-func-row-desc')?.textContent ?? ''};
});

async function checkPickerSearch(page) {
  await ensureDemog(page);
  await openSpec(page, SEARCH_SPEC);
  await waitStatus(page, 'searchForm');
  await expandTree(page);
  await selectRow(page, 'searchField');
  await openSection(page, 'Bindings');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="value"]').first().click();
  await page.waitForSelector('.u2-bind-picker', {timeout: 5000});
  await pickerSearch(page, '.u2-bind-picker').fill('subj');
  await page.waitForTimeout(400);
  const revealed = await pickerLabels(page);
  await shot(page, 'pickers-2-revealed');
  ok('pickers/2a/a-search-reveals-its-match-instead-of-leaving-it-under-closed-branches',
    revealed.some((label) => label.startsWith('subj')) &&
    revealed.some((label) => label.startsWith('currentRow')), JSON.stringify(revealed));

  await pickerSearch(page, '.u2-bind-picker').fill('');
  await page.waitForTimeout(400);
  const cleared = await pickerLabels(page);
  await shot(page, 'pickers-2-cleared');
  ok('pickers/2b/clearing-the-search-keeps-open-what-it-opened',
    cleared.some((label) => label.startsWith('subj')) &&
    cleared.some((label) => label.startsWith('df ')), JSON.stringify(cleared));

  await pickerSearch(page, '.u2-bind-picker').fill('u2Record');
  await page.waitForTimeout(400);
  const noBinding = await page.locator('.u2-bind-picker .u2-picker-empty').first().textContent();
  await shot(page, 'pickers-2-no-binding');
  ok('pickers/2c/an-empty-binding-result-says-so', /No bindings match "u2Record"/.test(noBinding ?? ''),
    `"${noBinding}"`);
  await dialogCancel(page, '.u2-bind-picker');
  await page.waitForTimeout(300);

  await selectRow(page, 'searchButton');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="click"]').first().click();
  await page.waitForSelector('.u2-func-picker', {timeout: 5000});
  await pickerSearch(page, '.u2-func-picker').fill('zzNoSuchFunction');
  await page.waitForTimeout(400);
  const noFunc = await page.locator('.u2-func-picker .u2-picker-empty').first().textContent();
  await shot(page, 'pickers-2-no-function');
  ok('pickers/2d/an-empty-function-result-says-so-instead-of-a-blank-white-box',
    /No functions match "zzNoSuchFunction"/.test(noFunc ?? ''), `"${noFunc}"`);

  await pickerSearch(page, '.u2-func-picker').fill('demoOrders');
  await page.waitForTimeout(400);
  const layout = await funcRowLayout(page);
  const names = await page.evaluate(() =>
    [...document.querySelectorAll('.u2-func-picker .u2-list-row')].map((row) => row.dataset.u2Func));
  await shot(page, 'pickers-2-func-row');
  ok('pickers/2e/the-row-lays-its-parts-out-and-shows-what-the-search-matched-on',
    layout !== null && layout.gap >= 0 && layout.top >= 0 && layout.bottom >= 0 &&
    layout.description !== '' && names.every((name) => /^[A-Za-z_$][\w$]*$/.test(name ?? '')),
    `${JSON.stringify(layout)} names=${JSON.stringify(names)}`);
  await dialogCancel(page, '.u2-func-picker');
  await page.waitForTimeout(300);
}

export const checks = [
  {id: 'pickers/1 binding picker: the Add-binding row, the tree walk, undo, the cycle refusal', run: checkBindPicker},
  {id: 'pickers/2 picker search: the match revealed, the expansion kept, the empty states, the row layout',
    run: checkPickerSearch},
];
