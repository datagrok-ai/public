/* The gestures steps are made of. Each one takes an element phrase and encodes the platform's
   quirks once: real pointer clicks for canvases, key-by-key typing (Dart change listeners ignore
   `fill`), native `<select>` first for choices, the editor part of a composite input. */
import {resolve} from 'node:path';
import {type Locator, type Page} from '@playwright/test';
import {expect} from './patience.js';
import type {ElementRef} from './args.js';
import {cssString, escapeRegExp, exactText, locateActionable as locate, refOf, withAttr} from './locate.js';

// a real control first: a Dart float input puts a `div.ui-input-editor` wrapper before its
// `input.ui-input-editor`, and `.first()` takes DOM order, so a bare `.ui-input-editor` in the same
// list would hand back the div — "Not an input element" on the first value read
const CONTROL = '[data-u2-part="editor"] input, [data-u2-part="editor"] select, [data-u2-part="editor"] textarea, ' +
  '[data-u2-part="editor"][contenteditable], input, select, textarea, [contenteditable="true"]';
const EDITOR = '.ui-input-editor';
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
  for (const selector of [CONTROL, EDITOR, EDITOR_PART]) {
    const inner = loc.locator(selector).first();
    if (await inner.count() > 0)
      return inner;
  }
  return loc;
}

export async function hasFocus(loc: Locator): Promise<boolean> {
  return loc.evaluate((e) => e === document.activeElement || e.contains(document.activeElement)).catch(() => false);
}

/** Types the text over whatever the editor holds, and again until it holds exactly that. A
 * keystroke that creates or rebuilds the editor lands at an unpredictable moment ("SEX" came back
 * "EXS"), and an editor the widget closes under the typing keeps only what came after it (a grid
 * cell typed "51" and committed "1"). An editor whose value cannot be read back is typed into once,
 * and so is one that refuses the text (read-only, disabled): refusing it is what the step after
 * such a typing claims. */
export async function typeVerified(editor: Locator, text: string, what: string): Promise<void> {
  await editor.press('Control+A');
  await editor.pressSequentially(text);
  const refuses = await editor.evaluate((e) => (e as HTMLInputElement).readOnly || (e as HTMLInputElement).disabled ||
    e.getAttribute('aria-readonly') === 'true' || e.getAttribute('aria-disabled') === 'true').catch(() => false);
  if (refuses || await editor.inputValue().catch(() => null) === null)
    return;
  await expect.poll(async () => {
    if (await editor.inputValue() !== text) {
      await editor.press('Control+A');
      await editor.pressSequentially(text);
    }
    return editor.inputValue();
  }, {timeout: 5000, message: `the text typed into ${what}`}).toBe(text);
}

export async function typeInto(page: Page, target: ElementRef, text: string, commit = false): Promise<void> {
  const editor = await editorOf(page, target);
  // an editor that has the focus already (a cell editor the double-click just opened) is not
  // clicked: a click is what would blur it, and a blurred cell editor commits and closes
  if (!await hasFocus(editor))
    await editor.click();
  if (gestureOf(page, target).type === 'fill')
    await editor.fill(text);
  else
    await typeVerified(editor, text, target.phrase);
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

/** Types a column name into the open picker and presses Enter, without claiming the pick landed —
 * a selector that does not offer that column keeps the one it had, and the feature reads it
 * afterwards. Returns the popup, for the caller that does claim it. */
export async function typeInColumnGrid(page: Page, option: string, what: string, selector?: Locator): Promise<Locator> {
  const popup = page.locator('.d4-column-grid').last();
  await popup.waitFor({state: 'visible', timeout: 10000});
  // pressed ON the selector where there is one: it does not always keep the focus its own
  // mouse-down gave it, and a letter that lands elsewhere opens no search box at all
  await (selector ? selector.press(option[0]) : page.keyboard.press(option[0]));
  const search = page.locator('input.d4-column-selector-search-input');
  await search.waitFor({state: 'visible', timeout: 10000});
  // and the letter lands in it at an unpredictable moment — before anything typed after it, or
  // after all of it ("SEX" came back "EXS", "RACE" as "RACER") — so the name goes in over
  // whatever is there, until the box holds it and nothing else
  await typeVerified(search, option, `the column picker of ${what}`);
  await search.press('Enter');
  return popup;
}

/** Types a column name into the picker a `.d4-column-grid` popup opens, and commits it. The caller
 * opens the popup and decides where the pointer may go (the filter panel's picker is part of a
 * header that is only shown while the panel is hovered).
 *
 * Three things about this picker cost a run each to learn. Its search box does not exist until a
 * letter is typed at the selector, and the letter that creates it lands in it out of order with
 * anything typed while it was being created — "SEX" came back "EXS" — so the name is retyped over
 * whatever landed. Enter is pressed ON the box, because the grid moves the focus while it filters
 * and a keystroke sent to the page then goes nowhere. And a picker still open after Enter took no
 * column: it is left showing the row the pointer previewed, which is what a later assertion would
 * otherwise report as the wrong column. */
export async function pickInColumnGrid(page: Page, option: string, what: string, selector?: Locator): Promise<void> {
  const popup = await typeInColumnGrid(page, option, what, selector);
  await popup.waitFor({state: 'detached', timeout: 5000}).catch(() => {
    throw new Error(`the column picker of ${what} is still open after Enter: "${option}" was not taken`);
  });
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
    // the popup opens under the pointer, and the row it rests on is previewed onto the selector,
    // which is what a picker that took nothing is left showing
    await page.mouse.move(2, 2);
    await pickInColumnGrid(page, option, target.phrase, columnSelector);
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
  // the Dart switch keeps its checkbox hidden and takes the click on a div, so it is not settable
  if (await loc.first().locator('.ui-input-switch').count() > 0)
    return setSwitched(page, target, checked);
  const box = loc.locator('input[type="checkbox"], input[type="radio"], [role="checkbox"], [role="switch"]').first();
  if (await box.count() > 0) {
    await box.setChecked(checked);
    return;
  }
  await loc.setChecked(checked);
}

export async function toggle(page: Page, target: ElementRef): Promise<void> {
  const loc = await locate(page, target);
  const box = loc.locator('.ui-input-switch, input[type="checkbox"], [role="checkbox"], [role="switch"]').first();
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
/** Where an element says it is open: `aria-expanded` on itself or on its header inside, and — for
 * the Dart tree, which has neither — the class its twistie carries. `null` when it says nothing. */
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

/** Reading the state first is what makes this idempotent: without it, a tree row that is already
 * open is closed by "user expands". */
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

/** A switch that governs something: the Dart `SwitchInput` draws it as `div.ui-input-switch` and
 * keeps its real checkbox hidden, so neither a click nor a read can go through the control. */
const SWITCH = '[role="switch"], .ui-input-switch';

/** The switch of an element, which a parameter form does not put inside it: the sensitivity
 * analysis form inserts the switch into the input's own host, the fitting form leaves it as a
 * separate input before it (`getSwitchElement`, compute-utils). A part selector cannot climb out
 * of the element and a composition would name the wrong parameter's switch, so the search is one
 * gesture — itself, inside, then back over the siblings — and the one it finds is marked, the way
 * a typed-into editor is, because a sibling has no phrase of its own. */
export async function switchOf(page: Page, target: ElementRef): Promise<Locator> {
  const self = (await locate(page, target)).first();
  const found = await self.evaluate((el, sel) => {
    for (const old of Array.from(document.querySelectorAll('[data-bdd-switch]')))
      old.removeAttribute('data-bdd-switch');
    let hit: Element | null = el.matches(sel) ? el : el.querySelector(sel);
    for (let sib = el.previousElementSibling; hit === null && sib !== null; sib = sib.previousElementSibling)
      hit = sib.matches(sel) ? sib : sib.querySelector(sel);
    if (hit === null)
      return false;
    hit.setAttribute('data-bdd-switch', '');
    return true;
  }, SWITCH);
  if (!found)
    throw new Error(`${target.phrase} has no switch, inside it or beside it`);
  return page.locator('[data-bdd-switch]');
}

/** `aria-checked` when the platform says it, else the class the Dart switch carries. */
export function readSwitch(loc: Locator): Promise<boolean | null> {
  return loc.first().evaluate((el) => {
    const aria = el.getAttribute('aria-checked');
    if (aria !== null)
      return aria === 'true';
    const sw = el.classList.contains('ui-input-switch') ? el : el.querySelector('.ui-input-switch');
    return sw === null ? null : sw.classList.contains('ui-input-switch-on');
  });
}

/** Idempotent, like `setExpanded`: a switch already on stays on. */
export async function setSwitched(page: Page, target: ElementRef, on: boolean): Promise<void> {
  const sw = await switchOf(page, target);
  if (await readSwitch(sw) === on)
    return;
  await sw.click();
  await expect.poll(() => readSwitch(sw),
    {message: `the switch of ${target.phrase} after switching it ${on ? 'on' : 'off'}`}).toBe(on);
}

/** A line at the top of a code editor. An editor's document is not an input value — CodeMirror
 * keeps it in its own model behind a hidden textarea — so the text goes in through the keyboard at
 * the caret, and the claim is the text the editor then shows. Control+Home rather than a click at
 * the first line: an editor scrolled down renders no first line to click. */
export async function insertLine(page: Page, target: ElementRef, text: string): Promise<void> {
  const loc = (await locate(page, target)).first();
  await loc.click();
  await page.keyboard.press('Control+Home');
  await page.keyboard.type(text);
  await page.keyboard.press('Enter');
  await expect(loc, `${target.phrase} after the line was typed`).toContainText(text);
}

/** A value set by dragging the slider of an input, not by typing into it: a real pointer press on
 * the thumb, a walk to where the value lives on the track, and a release. The track says what it
 * spans (`min`, `max`, `step`), so a feature can name the value a reader would aim at; a pixel of a
 * 200-px track is worth a hundredth of the range, so the landing is checked against the coarser of
 * one step and one pixel. */
export async function dragSlider(page: Page, target: ElementRef, value: number): Promise<void> {
  const loc = (await locate(page, target)).first();
  const range = loc.locator('input[type="range"]').first();
  if (await range.count() === 0)
    throw new Error(`${target.phrase} has no slider to drag`);
  const track = await range.evaluate((el: HTMLInputElement) =>
    ({min: Number(el.min), max: Number(el.max), step: Number(el.step) || 0, value: Number(el.value)}));
  if (value < track.min || value > track.max)
    throw new Error(`the slider of ${target.phrase} spans ${track.min} to ${track.max}, so it cannot be dragged to ${value}`);
  const box = await range.boundingBox();
  if (box === null)
    throw new Error(`the slider of ${target.phrase} has no rectangle on the page`);
  const at = (v: number): number => box.x + box.width * ((v - track.min) / (track.max - track.min));
  const y = box.y + box.height / 2;
  const from = at(track.value);
  const to = at(value);
  await page.mouse.move(from, y);
  await page.mouse.down();
  for (let i = 1; i <= 8; i++)
    await page.mouse.move(from + (to - from) * i / 8, y);
  await page.mouse.up();
  const tolerance = Math.max(track.step, (track.max - track.min) / box.width) * 2;
  await expect.poll(async () => Math.abs(Number(await range.inputValue()) - value) <= tolerance,
    {message: `the slider of ${target.phrase} within ${tolerance} of ${value} after the drag`}).toBe(true);
}
