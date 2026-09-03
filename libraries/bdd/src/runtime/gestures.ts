/* The gestures steps are made of. Each one takes an element phrase and encodes the platform's
   quirks once: real pointer clicks for canvases, key-by-key typing (Dart change listeners ignore
   `fill`), native `<select>` first for choices, the editor part of a composite input. */
import type {Locator, Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {cssString, escapeRegExp, exactText, locate, refOf, withAttr} from './locate.js';

const EDITOR = '[data-u2-part="editor"] input, [data-u2-part="editor"] select, [data-u2-part="editor"] textarea, ' +
  '[data-u2-part="editor"][contenteditable], .ui-input-editor, input, select, textarea, [contenteditable="true"]';
// popup triggers (icon, function and columns pickers) are the editor part itself
const EDITOR_PART = '[data-u2-part="editor"]';
const OPTION = '[role="option"], .u2-menu-item, .u2-combobox-option, .d4-menu-item, .u2-list-row';
const OPTION_LABEL = '.u2-fb-label, .u2-typeahead-text, .u2-typeahead-user-name, .u2-multi-select-text, .u2-menu-label';
const CLOSE = '.u2-dialog-close, [aria-label="Close"], [name="icon-times"], .d4-dialog-close';
const TWISTIE = '.u2-tree-twistie, .d4-tree-view-tri';

function gestureOf(page: Page, target: ElementRef): {click?: 'mouse' | 'dom'; type?: 'keyboard' | 'fill'} {
  const ref = refOf(page, target);
  const plan = ref.plan;
  return (plan.type === 'entry' ? plan.entry.gestures : plan.type === 'kind' ? plan.kind.gestures : undefined) ?? {};
}

export async function click(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  if (gestureOf(page, target).click === 'mouse') {
    await loc.scrollIntoViewIfNeeded();
    const box = await loc.boundingBox();
    if (!box)
      throw new Error(`${target.phrase}: no bounding box to click`);
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
    return;
  }
  await loc.click();
}

export async function dblclick(page: Page, target: ElementRef): Promise<void> {
  await (await locate(page, target)).dblclick();
}

export async function rightclick(page: Page, target: ElementRef): Promise<void> {
  await (await locate(page, target)).click({button: 'right'});
}

/** Leaves the element first — a pointer already resting on it (the previous click) produces no
 * pointerenter, and tooltips listen for that — then stays until the layout has settled: a shift
 * under the pointer right after the move (a view still docking) leaves it again, unseen. */
export async function hover(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  await loc.scrollIntoViewIfNeeded();
  for (let attempt = 0; attempt < 3; attempt++) {
    const before = await loc.boundingBox();
    if (before)
      await page.mouse.move(Math.max(0, before.x - 8), Math.max(0, before.y - 8));
    await loc.hover();
    await page.waitForTimeout(350);
    const after = await loc.boundingBox();
    if (!before || !after || (before.x === after.x && before.y === after.y))
      return;
  }
}

/** The editable control of an element: itself when it is one, otherwise its editor part. */
export async function editorOf(page: Page, target: ElementRef): Promise<Locator> {
  const loc = await locate(page, target);
  const editable = await loc.evaluate((e) => ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.tagName) ||
    (e as HTMLElement).isContentEditable).catch(() => false);
  if (editable)
    return loc;
  for (const selector of [EDITOR, EDITOR_PART]) {
    const inner = loc.locator(selector).first();
    if (await inner.count() > 0)
      return inner;
  }
  return loc;
}

export async function typeInto(page: Page, target: ElementRef, text: string, commit = false): Promise<void> {
  const editor = await editorOf(page, target);
  await editor.click();
  if (gestureOf(page, target).type === 'fill') {
    await editor.fill(text);
  }
  else {
    await editor.press('Control+A');
    await editor.pressSequentially(text);
  }
  if (commit)
    await editor.press('Tab');
}

export async function clear(page: Page, target: ElementRef): Promise<void> {
  const editor = await editorOf(page, target);
  await editor.click();
  await editor.press('Control+A');
  await editor.press('Delete');
}

export function normalizeKey(key: string): string {
  const names: Record<string, string> = {ctrl: 'Control', control: 'Control', cmd: 'Meta', meta: 'Meta',
    alt: 'Alt', shift: 'Shift', esc: 'Escape', escape: 'Escape', enter: 'Enter', return: 'Enter',
    tab: 'Tab', space: 'Space', backspace: 'Backspace', delete: 'Delete', del: 'Delete',
    up: 'ArrowUp', down: 'ArrowDown', left: 'ArrowLeft', right: 'ArrowRight', home: 'Home', end: 'End',
    pageup: 'PageUp', pagedown: 'PageDown'};
  return key.split('+').map((part) => {
    const p = part.trim();
    return names[p.toLowerCase()] ?? (p.length === 1 ? p.toUpperCase() : p);
  }).join('+');
}

export async function press(page: Page, key: string): Promise<void> {
  await page.keyboard.press(normalizeKey(key));
}

export async function select(page: Page, target: ElementRef, option: string): Promise<void> {
  const loc = await locate(page, target);
  const native = loc.locator('select').first();
  if (await native.count() > 0) {
    await native.selectOption({label: option});
    return;
  }
  const editor = await editorOf(page, target);
  await editor.click();
  let options = optionsNamed(page, option);
  await options.first().waitFor({timeout: 1500}).catch(() => undefined);
  // comboboxes and typeaheads open on a keystroke, not on the click
  if (await options.count() === 0 && await editor.getAttribute('role') === 'combobox') {
    await editor.press('ArrowDown');
    await options.first().waitFor({timeout: 1500}).catch(() => undefined);
  }
  if (await options.count() === 0)
    options = page.locator(OPTION).filter({hasText: new RegExp(escapeRegExp(option), 'i')});
  await options.first().click();
}

/** Popup rows called `option`: by their whole text, by their primary-text part, or by a title/aria label. */
function optionsNamed(page: Page, option: string): Locator {
  const all = page.locator(OPTION);
  const exact = exactText(option);
  return all.filter({hasText: exact})
    .or(all.filter({has: page.locator(OPTION_LABEL, {hasText: exact})}))
    .or(page.locator(withAttr(OPTION, `[title="${cssString(option)}" i]`)))
    .or(page.locator(withAttr(OPTION, `[aria-label="${cssString(option)}" i]`)));
}

export async function setChecked(page: Page, target: ElementRef, checked: boolean): Promise<void> {
  const loc = await locate(page, target);
  const box = loc.locator('input[type="checkbox"], input[type="radio"], [role="checkbox"], [role="switch"]').first();
  if (await box.count() > 0) {
    await box.setChecked(checked);
    return;
  }
  await loc.setChecked(checked);
}

export async function toggle(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  const box = loc.locator('input[type="checkbox"], [role="checkbox"], [role="switch"]').first();
  await (await box.count() > 0 ? box : loc).click();
}

export async function close(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  const button = loc.locator(CLOSE).first();
  if (await button.count() > 0)
    await button.click();
  else
    await page.keyboard.press('Escape');
}

export async function focus(page: Page, target: ElementRef): Promise<void> {
  await (await editorOf(page, target)).focus();
}

export async function pressIn(page: Page, target: ElementRef, key: string): Promise<void> {
  await (await editorOf(page, target)).press(normalizeKey(key));
}

export async function drag(page: Page, source: ElementRef, target: ElementRef): Promise<void> {
  await (await locate(page, source)).dragTo(await locate(page, target));
}

export async function scrollTo(page: Page, target: ElementRef): Promise<void> {
  await (await locate(page, target)).first().scrollIntoViewIfNeeded();
}

/** Trees, accordion panes and dropdowns say where they are through `aria-expanded` — on the
 * element itself or on its header/trigger inside. */
export async function setExpanded(page: Page, target: ElementRef, expanded: boolean): Promise<void> {
  const self = (await locate(page, target)).first();
  const inner = self.locator('[aria-expanded]').first();
  const control = await self.getAttribute('aria-expanded') !== null ? self : await inner.count() > 0 ? inner : self;
  if (await control.getAttribute('aria-expanded') === String(expanded))
    return;
  // a tree row selects on click and toggles on its twistie
  const twistie = control.locator(TWISTIE).first();
  await (await twistie.count() > 0 ? twistie : control).click();
}
