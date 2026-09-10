/* The gestures steps are made of. Each one takes an element phrase and encodes the platform's
   quirks once: real pointer clicks for canvases, key-by-key typing (Dart change listeners ignore
   `fill`), native `<select>` first for choices, the editor part of a composite input. */
import {resolve} from 'node:path';
import {expect, type Locator, type Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {cssString, escapeRegExp, exactText, locateActionable as locate, refOf, withAttr} from './locate.js';

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

/** A click with a key or a chord held (Control adds to a selection, Shift extends it). */
export async function clickHolding(page: Page, target: ElementRef, key: string): Promise<void> {
  const keys = normalizeKey(key).split('+');
  for (const k of keys)
    await page.keyboard.down(k);
  try {
    await click(page, target);
  }
  finally {
    for (const k of keys.reverse())
      await page.keyboard.up(k);
  }
}

/** Puts the text on the page's clipboard and pastes it into the element's editor with Ctrl+V — the
 * platform's own paste handling runs (a comma list becomes an alternation, a newline list a union). */
export async function paste(page: Page, target: ElementRef, text: string): Promise<void> {
  const editor = await editorOf(page, target);
  await editor.click();
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  await page.evaluate((t) => navigator.clipboard.writeText(t), text);
  await page.keyboard.press('Control+V');
}

export async function dblclick(page: Page, target: ElementRef): Promise<void> {
  await (await locate(page, target)).dblclick();
}

export async function rightclick(page: Page, target: ElementRef): Promise<void> {
  await (await locate(page, target)).click({button: 'right'});
}

/** Clicks the element and answers the file chooser it opens with a file of the bdd project
 * (`file` relative to the project root, `BDD_ROOT`). */
export async function chooseFile(page: Page, target: ElementRef, file: string): Promise<void> {
  const loc = await locate(page, target);
  const chooser = page.waitForEvent('filechooser', {timeout: 5000});
  await loc.click();
  await (await chooser).setFiles(resolve(process.env.BDD_ROOT ?? process.cwd(), file));
}

/** The page's clipboard text — headless Chromium keeps a clipboard of its own per browser. */
export async function readClipboard(page: Page): Promise<string> {
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  return page.evaluate(() => navigator.clipboard.readText());
}

/** Leaves the element first, to its left on the same line — a pointer already resting on it (the
 * previous click) produces no pointerenter, and tooltips listen for that; leaving upwards would
 * cross a neighbouring menu row and close the submenu the element sits in — then lands on its
 * centre in one move (every pointer event costs a frame, and a Dart menu group opens on the first
 * move since 2026-09-07) and checks that the element is still where it was: a shift under the
 * pointer right after the move (a view still docking) leaves it again, unseen. The browser
 * coalesces mouse moves queued while its main thread is busy, so the pair can collapse into the
 * last move alone — one that enters nothing when the pointer already rested inside; the gesture
 * therefore waits, in the page and for a few frames at most, for the element's own `mouseenter`,
 * and repeats the pair when it did not come. */
export async function hover(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  const viewport = page.viewportSize();
  for (let attempt = 0; attempt < 3; attempt++) {
    let before = await loc.boundingBox();
    if (before && viewport && (before.y < 0 || before.x < 0 || before.y + before.height > viewport.height || before.x + before.width > viewport.width)) {
      await loc.scrollIntoViewIfNeeded();
      before = await loc.boundingBox();
    }
    if (!before) {
      await loc.hover();
      return;
    }
    const cy = before.y + before.height / 2;
    await loc.evaluate((el) => {
      (el as any).__bddEntered = false;
      el.addEventListener('mouseenter', () => { (el as any).__bddEntered = true; }, {once: true});
    });
    await page.mouse.move(Math.max(0, before.x - 8), cy);
    await page.mouse.move(before.x + before.width / 2, cy);
    const entered = await loc.evaluate((el) => new Promise<boolean>((resolve) => {
      let frames = 0;
      const tick = () => (el as any).__bddEntered || frames++ > 4 ? resolve((el as any).__bddEntered === true) : requestAnimationFrame(tick);
      tick();
    }));
    const after = await loc.boundingBox();
    if (entered && (!after || (before.x === after.x && before.y === after.y)))
      return;
  }
}

/** The editable control of an element: itself when it is one, otherwise its editor part. */
export async function editorOf(page: Page, target: ElementRef): Promise<Locator> {
  const loc = await locate(page, target);
  // a viewer handles keys on its root (a form walks rows on the arrows, a plot zooms on +/-);
  // its first input is a field, not its editor
  const own = await loc.evaluate((e) => ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.tagName) ||
    (e as HTMLElement).isContentEditable || e.matches('[name^="viewer-"], .d4-viewer')).catch(() => false);
  if (own)
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
  // the Dart column selector: a mouse-down opens its column grid, typing opens the grid's
  // search box, and Enter there takes the name typed as the column
  const columnSelector = (await loc.evaluate((el) => el.classList.contains('d4-column-selector'))) ? loc : loc.locator('.d4-column-selector').first();
  if (await columnSelector.count() > 0) {
    const box = await columnSelector.boundingBox();
    if (!box)
      throw new Error(`${target.phrase} has no box`);
    await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
    await page.mouse.down();
    await page.mouse.up();
    await page.locator('.d4-column-grid').last().waitFor({state: 'visible', timeout: 5000});
    await page.keyboard.type(option);
    await page.keyboard.press('Enter');
    await expect(columnSelector.locator('.d4-column-selector-column')).toHaveText(exactText(option), {timeout: 5000});
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

/** Where an element says it is open: `aria-expanded` on itself or on its header/trigger inside,
 * and — for the Dart tree, which has neither — the class its twistie carries. `null` when the
 * element says nothing (a leaf row has no twistie), so a caller can report that instead of
 * guessing. */
export function readExpanded(loc: Locator): Promise<boolean | null> {
  return loc.first().evaluate((el) => {
    const aria = el.getAttribute('aria-expanded') ?? el.querySelector('[aria-expanded]')?.getAttribute('aria-expanded');
    if (aria != null)
      return aria === 'true';
    const twistie = el.matches('.d4-tree-view-tri, .u2-tree-twistie') ? el :
      el.querySelector('.d4-tree-view-tri, .u2-tree-twistie');
    return twistie === null ? null :
      twistie.classList.contains('d4-tree-view-tri-expanded') || twistie.classList.contains('u2-tree-twistie-expanded');
  });
}

/** Trees, accordion panes and dropdowns say where they are through `aria-expanded` — on the
 * element itself or on its header/trigger inside — and the Dart tree through its twistie's class.
 * Reading it first is what makes this idempotent: without it, a tree row that is already open is
 * closed by "user expands". */
export async function setExpanded(page: Page, target: ElementRef, expanded: boolean): Promise<void> {
  const self = (await locate(page, target)).first();
  if (await readExpanded(self) === expanded)
    return;
  const inner = self.locator('[aria-expanded]').first();
  const control = await self.getAttribute('aria-expanded') !== null ? self : await inner.count() > 0 ? inner : self;
  // a tree row selects on click and toggles on its twistie
  const twistie = control.locator(TWISTIE).first();
  await (await twistie.count() > 0 ? twistie : control).click();
  await expect.poll(() => readExpanded(self),
    {message: `${target.phrase} after ${expanded ? 'expanding' : 'collapsing'} it`}).toBe(expanded);
}
