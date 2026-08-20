/* P3.5 WO-5 — multi-selection (F5): Ctrl/Shift gestures on the tree, Ctrl on the canvas, the
   multi panel, the compound delete with its single undo entry, and Esc collapsing to the lead.
   Then the second-acceptance fixes: a drop on the canvas padding lands in the root, a dblclick
   anywhere on a fresh row renames, and the guarded Events field commits on change. */
import {ok, shot} from '../local.mjs';
import {CAPTION, EDIT_SPEC, balloons, caption, clearBalloons, confirmDiscard, dumpViaDialog, expandTree,
  named, openSpec, panelField, ribbon, row, selectRow, statusText, surfaceNode, waitCaption, waitStatus}
  from '../lib.mjs';

export async function fixture(page) {
  await page.locator('.u2-palette input').first().fill('');
  await page.waitForTimeout(150);
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
}

const visibleAdorners = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.u2-designer-selected')]
    .filter((box) => getComputedStyle(box).display === 'block').length);

/** The multi status line is `'N selected · M nodes'` — no node path to wait on by name. */
async function waitSelected(page, n) {
  await page.waitForFunction((n) => [...document.querySelectorAll('div')]
    .some((x) => x.children.length === 0 && new RegExp(`^${n} selected · \\d+ node`).test(x.textContent)),
  n, {timeout: 3000}).catch(() => {});
  return statusText(page);
}

async function checkMultiSelect(page) {
  await selectRow(page, 'nameInput');
  // the platform panel coalesces shell.o assignments below ~1-2 s apart (the known platform
  // defect, view-and-panel/7d) — the multi gestures are spaced past that window so each one renders
  await page.waitForTimeout(2000);
  await row(page, 'saveButton').click({modifiers: ['Control']});
  const two = await waitSelected(page, 2);
  const waitNodes = (text) => page.waitForFunction(({sel, text}) =>
    (document.querySelector(sel)?.textContent ?? '') === text,
  {sel: CAPTION, text}, {timeout: 4000}).then(() => true, () => false);
  let converged = await waitNodes('2 nodes');
  if (!converged) {
    // toggle off and back on, each its own well-spaced user gesture — same two members after
    await page.waitForTimeout(2000);
    await row(page, 'saveButton').click({modifiers: ['Control']});
    await page.waitForTimeout(2000);
    await row(page, 'saveButton').click({modifiers: ['Control']});
    converged = await waitNodes('2 nodes');
  }
  await shot(page, 'multi-select-and-guards-1-multi-select');
  ok('multi-select-and-guards/1a/ctrl-click-selects-two', /^2 selected · 7 nodes/.test(two) &&
    await visibleAdorners(page) === 2 && converged,
  `status="${two}" adorners=${await visibleAdorners(page)} caption="${await caption(page)}"`);

  // the anchor is the last plain click — nameInput — so the range replaces the set
  await row(page, 'tabs').click({modifiers: ['Shift']});
  const three = await waitSelected(page, 3);
  ok('multi-select-and-guards/1b/shift-click-ranges-from-the-anchor', /^3 selected · 7 nodes/.test(three) &&
    await visibleAdorners(page) === 3, `status="${three}" adorners=${await visibleAdorners(page)}`);

  const box = await surfaceNode(page, 'heading').boundingBox();
  await page.keyboard.down('Control');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.keyboard.up('Control');
  const four = await waitSelected(page, 4);
  ok('multi-select-and-guards/1c/canvas-ctrl-click-adds-a-member', /^4 selected · 7 nodes/.test(four) &&
    await visibleAdorners(page) === 4, `status="${four}" adorners=${await visibleAdorners(page)}`);

  // the cover: nameInput and firstPane ride the tabs remove; heading is the lead
  await page.keyboard.press('Delete');
  await page.waitForTimeout(300);
  const afterDelete = await statusText(page);
  const gone = await named(page);
  await shot(page, 'multi-select-and-guards-1-multi-delete');
  ok('multi-select-and-guards/1d/multi-delete-removes-the-cover-and-selects-the-parent',
    ['tabs', 'nameInput', 'heading'].every((n) => !gone.includes(n)) && gone.includes('saveButton') &&
    /^editRoot · 2 nodes · modified/.test(afterDelete),
    `status="${afterDelete}" named=${JSON.stringify(gone)}`);

  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  const restored = await named(page);
  const undone = await statusText(page);
  ok('multi-select-and-guards/1e/one-undo-restores-the-whole-compound',
    ['heading', 'tabs', 'nameInput', 'saveButton'].every((n) => restored.includes(n)) &&
    !/modified/.test(undone) && /7 nodes/.test(undone),
    `status="${undone}" named=${JSON.stringify(restored)}`);

  await selectRow(page, 'nameInput');
  await row(page, 'saveButton').click({modifiers: ['Control']});
  await waitSelected(page, 2);
  await page.keyboard.press('Escape');
  const collapsed = await waitStatus(page, 'saveButton');
  ok('multi-select-and-guards/1f/esc-collapses-to-the-lead', /› saveButton · 7 nodes/.test(collapsed) &&
    !/selected/.test(collapsed) && await visibleAdorners(page) === 1,
  `status="${collapsed}" adorners=${await visibleAdorners(page)}`);

  // the 2b regression, with the multi machinery in place: a plain re-click still re-asserts —
  // converged with the same click-again pattern every single-select check uses (selectRow)
  await page.evaluate(() => {
    grok.shell.o = grok.shell.user;
  });
  await page.waitForTimeout(1200);
  await row(page, 'saveButton').click();
  let back = await waitCaption(page, 'saveButton', 2500);
  if (!back.startsWith('saveButton (')) {
    await page.waitForTimeout(1200);
    await row(page, 'saveButton').click();
    back = await waitCaption(page, 'saveButton');
  }
  await shot(page, 'multi-select-and-guards-1-reclick');
  ok('multi-select-and-guards/1g/plain-reclick-on-the-lead-reasserts', back === 'saveButton (u2-button)',
    `caption="${back}"`);
}

/** A refused value warns with the complete message and keeps the typed text. Builds its own
 * scratch form off the blank sample. */
async function checkPaddingDropAndGuardedFields(page) {
  await ribbon(page, 'New form').click();
  await page.waitForTimeout(300);
  if (await page.locator('.u2-dialog').count() > 0)
    await confirmDiscard(page);
  await waitStatus(page, 'form1');

  // the drop point is the surface padding well below the empty form — no rendered node there
  const surface = await page.locator('.u2-designer-surface').boundingBox();
  const form = await page.locator('.u2-designer-surface [data-u2-name="form1"]').boundingBox();
  const item = page.locator('.u2-palette-item[data-u2-tag="u2-text-input"]').first();
  const a = await item.boundingBox();
  await page.mouse.move(a.x + a.width / 2, a.y + a.height / 2);
  await page.mouse.down();
  const x = form.x + form.width / 2;
  const y = Math.min(form.y + form.height + 60, surface.y + surface.height - 5);
  await page.mouse.move(x, y, {steps: 25});
  await page.mouse.move(x, y + 1, {steps: 3});
  await page.mouse.up();
  await waitStatus(page, 'textInput1');
  const status = await statusText(page);
  await shot(page, 'multi-select-and-guards-2-padding-drop');
  ok('multi-select-and-guards/2a/canvas-padding-drop-inserts-into-the-root',
    /textInput1 · 2 nodes · modified/.test(status),
    `status="${status}" drop at ${Math.round(x)},${Math.round(y)} form bottom ${Math.round(form.y + form.height)}`);

  // Jane's gesture: a dblclick aimed at the ROW — its center sits right of the short label
  await row(page, 'textInput1').dblclick();
  const editor = page.locator('.u2-tree-rename');
  const renaming = await editor.count() === 1 && await editor.inputValue() === 'textInput1';
  ok('multi-select-and-guards/2b/row-dblclick-opens-rename-past-the-label', renaming, `editor=${await editor.count()}`);
  await page.keyboard.press('Escape');

  await selectRow(page, 'textInput1');
  const eventField = panelField(page, 'input');
  await eventField.click();
  await clearBalloons(page);
  await eventField.pressSequentially('cmd', {delay: 30});
  await page.waitForTimeout(300);
  const midBalloons = await balloons(page);
  const midText = await eventField.inputValue();
  ok('multi-select-and-guards/2c/partial-typing-neither-commits-nor-wipes',
    midBalloons.length === 0 && midText === 'cmd',
    `balloons=${JSON.stringify(midBalloons)} field="${midText}"`);

  await eventField.pressSequentially(':ping', {delay: 30});
  await eventField.press('Enter');
  await page.waitForTimeout(300);
  const dump = await dumpViaDialog(page);
  ok('multi-select-and-guards/2d/the-full-value-commits-on-enter',
    JSON.parse(dump).root.children[0].on?.input === 'cmd:ping' &&
    (await balloons(page)).length === 0, dump.slice(0, 200));

  await eventField.fill('go');
  await eventField.press('Enter');
  await page.waitForFunction(() => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => /an event must name a command/.test(b.textContent)), null, {timeout: 5000}).catch(() => {});
  const refused = (await balloons(page)).find((b) => /an event must name a command/.test(b)) ?? '';
  await shot(page, 'multi-select-and-guards-2-refusal-keeps-text');
  ok('multi-select-and-guards/2e/a-refusal-warns-complete-and-keeps-the-text',
    refused.includes('an event must name a command as \'cmd:\' followed by the command name — got \'go\'') &&
    await eventField.inputValue() === 'go',
  `balloon="${refused}" field="${await eventField.inputValue()}"`);
  await clearBalloons(page);
}

export const checks = [
  {id: 'multi-select-and-guards/1 multi-selection: Ctrl/Shift, compound delete, Esc collapse', run: checkMultiSelect},
  {id: 'multi-select-and-guards/2 padding drops, row-dblclick rename, guarded event fields',
    run: checkPaddingDropAndGuardedFields},
];
