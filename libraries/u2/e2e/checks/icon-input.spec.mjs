/* The icon picker in the context panel: a tab pane's `icon` child prop and a button's `icon` prop
   render as IconInput, a pick lands as a set-prop patch the canvas and the dump both show, and
   one Ctrl+Z takes it back. */
import {ok, shot} from '../local.mjs';
import {EDIT_SPEC, dumpViaDialog, expandTree, focusTree, openSpec, panel, selectRow, waitStatus}
  from '../lib.mjs';

export async function fixture(page) {
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
}

const picker = (page, name) => page.locator(`.u2-designer-properties [data-u2-prop="${name}"] .u2-icon-input`).first();
const shownName = (page, name) => picker(page, name).locator('.u2-icon-input-name').textContent();

/** Opens the picker of one panel field, filters, and clicks the cell. */
async function pick(page, prop, icon) {
  await picker(page, prop).click();
  const popup = page.locator('.u2-icon-input-popup');
  await popup.waitFor({timeout: 5000});
  await popup.locator('input').fill(icon);
  await page.waitForTimeout(150);
  await popup.locator(`.u2-grid-cell[title="${icon}"]`).first().click();
  await page.waitForTimeout(400);
}

async function checkPaneIcon(page) {
  await selectRow(page, 'firstPane');
  const before = await panel(page);
  const field = await picker(page, 'icon').count();
  const empty = field === 1 ? await shownName(page, 'icon') : '';
  ok('icon-input/1a/child-icon-prop-renders-the-picker', field === 1 && 'icon' in (before['Parent (u2-tabs)'] ?? {}) && empty === '—',
    `pickers=${field} shown="${empty}" section=${JSON.stringify(before['Parent (u2-tabs)'] ?? null)}`);

  await pick(page, 'icon', 'table');
  await shot(page, 'icon-input-1-pane-icon');
  const shown = await shownName(page, 'icon');
  const glyph = await page.locator('.u2-designer-surface .u2-tabs-icon .fa-table').count();
  const dump = JSON.parse(await dumpViaDialog(page));
  const pane = dump.root.children[1].children[0];
  ok('icon-input/1b/pick-patches-the-child-prop', shown === 'table' && glyph === 1 && pane.props?.icon === 'table',
    `panel="${shown}" glyphs=${glyph} props=${JSON.stringify(pane.props)}`);

  await focusTree(page, 'firstPane');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(400);
  const undone = await page.locator('.u2-designer-surface .u2-tabs-icon .fa-table').count();
  ok('icon-input/1c/undo-takes-the-icon-back', undone === 0, `glyphs after undo=${undone}`);
}

async function checkButtonIcon(page) {
  await selectRow(page, 'saveButton');
  const field = await picker(page, 'icon').count();
  await pick(page, 'icon', 'save');
  await shot(page, 'icon-input-2-button-icon');
  const glyph = await page.locator('.u2-designer-surface [data-u2-name="saveButton"] .fa-save').count();
  const text = await page.locator('.u2-designer-surface [data-u2-name="saveButton"]').first().textContent();
  ok('icon-input/2a/button-icon-prop-renders-and-patches', field === 1 && glyph === 1 && text.trim() === 'Save',
    `pickers=${field} glyphs=${glyph} text="${text}"`);
}

export const checks = [
  {id: 'icon-input/1 tab pane icon through the picker', run: checkPaneIcon},
  {id: 'icon-input/2 button icon through the picker', run: checkButtonIcon},
];
