/* Outcome checks over Playwright's retrying `expect`, so a feature never needs a wait step. */
import {Locator, Page} from '@playwright/test';
import {expect} from './patience.js';
import type {ElementRef} from './args.js';
import {CHECKABLE, editorOf, openPropertyEditor, readExpanded, readSwitch, switchOf, whenPresent} from './gestures.js';
import {exactText, locate, locateActionable, refOf} from './locate.js';
import {DomReadings, noteContract, readControls, readState, reports, shownFirst, summary, textOf, U2State} from './u2-state.js';

import type {State} from '../states.js';

export type {State};
export {STATES} from '../states.js';
const INVALID_CLASSES = ['d4-invalid', 'd4-forced-invalid', 'u2-input-invalid'];

const ROWS = ['.u2-list-row', '[role="option"]', '[role="row"]', '[role="tab"]', 'option', '.d4-list-item', '[name="legend-item"]', 'tbody tr', 'tr', 'li'];
// a dock manager's tab, a gallery's view-mode icon and a Mol* viewport toggle say which is on with a
// class of their own; d4-current elsewhere in a gallery marks the current card, not a selection
const SELECTED = '[aria-selected="true"], [aria-pressed="true"], [aria-checked="true"], [aria-current]:not([aria-current="false"]), ' +
  '.u2-list-row-selected, .tab-handle-selected, .dockspan-tab-handle-selected, .grok-gallery-search-bar .d4-current, ' +
  '.msp-btn-link-toggle-on';

/** A state claim reads, per element, the status of a u2 control root that reports the state's key
 * (`u2-state.ts`), and today's DOM reading otherwise, over the same elements and quantifier as before
 * the status existed. */
export async function expectState(page: Page, target: ElementRef, state: State, negate = false): Promise<void> {
  const loc = ['visible', 'hidden', 'present', 'absent', 'enabled', 'disabled'].includes(state) ?
    await locate(page, target) : await locateActionable(page, target);
  const on: On = {page, phrase: target.phrase};
  switch (state) {
    case 'visible': return expectVisible(loc, !negate);
    case 'hidden': return expectVisible(loc, negate);
    case 'present': return negate ? expect(loc).toHaveCount(0) : expect(loc.first()).toBeAttached();
    case 'absent': return negate ? expect(loc.first()).toBeAttached() : expect(loc).toHaveCount(0);
    case 'enabled': return expectFlag(on, loc, 'disabled', negate);
    case 'disabled': return expectFlag(on, loc, 'disabled', !negate);
    case 'checked': return expectChecked(on, loc, !negate);
    case 'unchecked': return expectChecked(on, loc, negate);
    case 'partially checked': return expectMixed(on, loc, !negate);
    case 'invalid': return expectFlag(on, loc, 'invalid', !negate);
    case 'valid': return expectFlag(on, loc, 'invalid', negate);
    case 'ready': return expectFlag(on, loc, 'ready', !negate);
    case 'selected': return expectFlag(on, loc, 'selected', !negate);
    case 'empty': return expectEmpty(on, loc, target, !negate);
    case 'expanded': return expectExpanded(on, loc, !negate);
    case 'collapsed': return expectExpanded(on, loc, negate);
    case 'focused': return expectFocused(on, loc, !negate);
    case 'read-only': return expectFlag(on, loc, 'read-only', !negate);
    case 'editable': return expectFlag(on, loc, 'read-only', negate);
    case 'busy': return expectFlag(on, loc, 'busy', !negate);
    case 'required': return expectFlag(on, loc, 'required', !negate);
  }
}

interface On {
  page: Page;
  phrase: string;
}

type Flag = 'disabled' | 'invalid' | 'ready' | 'selected' | 'read-only' | 'busy' | 'required';

const FLAGS: Record<Flag, {words: [string, string]; attribute: string; message?: string}> = {
  'disabled': {words: ['disabled', 'enabled'], attribute: 'aria-disabled'},
  'invalid': {words: ['invalid', 'valid'], attribute: 'aria-invalid'},
  // a completed asynchronous widget explicitly clears aria-busy; missing markup is not readiness, and a
  // failed computation can stop running but must not count as a ready result
  'ready': {words: ['ready', 'not ready'], attribute: 'aria-busy', message: 'a completed, valid result (aria-busy=false)'},
  'selected': {words: ['selected', 'not selected'], attribute: 'aria-selected'},
  'read-only': {words: ['read-only', 'editable'], attribute: 'aria-readonly'},
  'busy': {words: ['busy', 'not busy'], attribute: 'aria-busy'},
  'required': {words: ['required', 'not required'], attribute: 'aria-required'},
};

/** In the page, per element: whether it is shown, today's DOM reading of the flag, and the status
 * reading when the element is a u2 control root whose status has the key. Disabled: the element or an
 * ancestor says so (`aria-disabled`, the u2/Dart disabled classes), or it — or the control inside it —
 * is natively disabled; a u2 control is also disabled inside a disabled ancestor. Invalid: `aria-invalid`
 * or the platform's invalid classes on the element or a control inside. Selected: options and tabs say
 * `aria-selected`, toggles and cards `aria-pressed`, radio-like buttons `aria-checked`, wizard steps and
 * breadcrumbs `aria-current`. Busy: the element or an ancestor. */
const flagReadings = (els: Element[], [flag, classes, selected]: [Flag, string[], string]) => {
  const U2 = (window as any).DG?.U2;
  return els.map((el) => {
    let s: any = null;
    try {
      s = U2?.control?.(el)?.root === el ? U2.state(el) : null;
    }
    catch {
      s = null;
    }
    const keys: string[] = s ? Object.keys(s) : [];
    const editor = (el.matches('input, select, textarea') ? el : (el.id ? document.getElementById(`${el.id}-editor`) : null) ??
      el.querySelector('input, select, textarea')) as HTMLInputElement | null;
    let dom = false;
    let status: boolean | undefined;
    switch (flag) {
      case 'disabled': {
        const chain: Element[] = [];
        for (let e: Element | null = el; e; e = e.parentElement)
          chain.push(e);
        const marked = chain.map((e) => e.getAttribute('aria-disabled') === 'true' ||
          ['u2-input-disabled', 'd4-disabled', 'd4-filter-disabled', 'd4-menu-item-disabled'].some((c) => e.classList.contains(c)));
        const above = marked.slice(1).some(Boolean);
        const native = el.matches('input, select, textarea, button, fieldset, option') ? el : el.querySelector('input, select, textarea, button');
        dom = marked.some(Boolean) || (native !== null && (native as HTMLInputElement).disabled === true);
        status = keys.includes('enabled') ? !s.enabled || above : undefined;
        break;
      }
      case 'invalid':
        dom = [el, ...Array.from(el.querySelectorAll('input, select, textarea, [aria-invalid]'))]
          .some((e) => e.getAttribute('aria-invalid') === 'true' || classes.some((c) => e.classList.contains(c)));
        status = keys.includes('valid') ? s.valid === false : undefined;
        break;
      case 'ready':
        dom = el.getAttribute('aria-busy') === 'false' && el.getAttribute('aria-invalid') !== 'true';
        status = keys.includes('busy') ? s.busy === false && s.valid !== false : undefined;
        break;
      case 'selected':
        dom = el.matches(selected);
        status = keys.includes('selected') || keys.includes('pressed') ? (s.selected ?? s.pressed) === true : undefined;
        break;
      case 'read-only':
        dom = el.getAttribute('aria-readonly') === 'true' || editor?.getAttribute('aria-readonly') === 'true' || editor?.readOnly === true;
        status = keys.includes('readOnly') ? s.readOnly === true : undefined;
        break;
      case 'busy':
        dom = el.closest('[aria-busy="true"]') !== null;
        status = keys.includes('busy') ? s.busy === true || el.parentElement?.closest('[aria-busy="true"]') != null : undefined;
        break;
      case 'required':
        dom = el.getAttribute('aria-required') === 'true' || editor?.getAttribute('aria-required') === 'true' || editor?.required === true;
        status = keys.includes('required') ? s.required === true : undefined;
        break;
    }
    return {shown: el.getClientRects().length > 0 && getComputedStyle(el).visibility !== 'hidden', dom, status};
  });
};

/** One query per poll. Disabled and the states added with the status hold over the visible matches when
 * there are any (the Dart menu's hidden mirror), else all of them (a property row in a panel that is not
 * shown still says whether it is gated); invalid and selected over any match; ready over every one. */
async function expectFlag(on: On, loc: Locator, flag: Flag, want: boolean): Promise<void> {
  const {words: [yes, no], attribute, message} = FLAGS[flag];
  const holds = async () => {
    const all = await loc.evaluateAll(flagReadings, [flag, INVALID_CLASSES, SELECTED] as [Flag, string[], string]);
    for (const r of all) {
      if (r.status !== undefined && r.status !== r.dom)
        noteContract(on.page, on.phrase, r.status ? yes : no, r.dom ? yes : no, attribute);
    }
    const value = (r: {dom: boolean; status?: boolean}) => r.status ?? r.dom;
    if (flag === 'invalid' || flag === 'selected')
      return all.some(value);
    if (flag === 'ready')
      return all.length > 0 && all.every(value);
    const shown = all.filter((r) => r.shown);
    const els = shown.length > 0 ? shown : all;
    return els.length === 0 ? undefined : els.every(value);
  };
  await expect.poll(holds, {message: message ?? `${want ? yes : no} expected`}).toBe(want);
}

/** The status claim of a state the DOM path checks through a Playwright matcher, on the first element:
 * polled, `undefined` while the element has no status (it is being rebuilt). */
async function expectStatus(on: On, loc: Locator, want: boolean, [yes, no]: [string, string], status: (s: U2State) => boolean,
  dom: (d: DomReadings) => boolean | null, attribute: string): Promise<void> {
  await expect.poll(async () => {
    const r = await readState(loc);
    if (!r?.state)
      return undefined;
    const s = status(r.state);
    const d = dom(r.dom);
    if (d !== null && d !== s)
      noteContract(on.page, on.phrase, s ? yes : no, d ? yes : no, attribute);
    return s;
  }, {message: `${want ? yes : no} expected`}).toBe(want);
}

/** The status says whether focus is inside the control or a popup it owns; the DOM path needs the element itself. */
async function expectFocused(on: On, loc: Locator, focused: boolean): Promise<void> {
  if (await reports(loc, ['focused'], true))
    return expectStatus(on, loc, focused, ['focused', 'not focused'], (s) => s.focused, () => null, 'focus');
  await (focused ? expect(loc) : expect(loc).not).toBeFocused();
}

/** The text the control shows (the value when it reports no text); the DOM path reads the editor. */
async function expectEmpty(on: On, loc: Locator, target: ElementRef, empty: boolean): Promise<void> {
  if (await reports(loc, ['text', 'value']))
    return expectStatus(on, loc, empty, ['empty', 'not empty'], (s) => typeof s.text === 'string' ? s.text === '' :
      s.value === null || (Array.isArray(s.value) && s.value.length === 0), (d) => d.value === null ? null : d.value === '', 'the editor\'s value');
  await (empty ? expect(await editorOf(on.page, target)) : expect(await editorOf(on.page, target)).not).toHaveValue('');
}

/** A parameter form's switch is not inside the input it governs, so this is a claim of its own
 * rather than a state of the element. */
export async function expectSwitched(page: Page, target: ElementRef, on: boolean): Promise<void> {
  const sw = await switchOf(page, target);
  await expect.poll(() => readSwitch(sw), {message: `the switch of ${target.phrase}`}).toBe(on);
}

/** Several matches (stacked balloons, repeated rows): visible when any is, hidden when none is —
 * one query either way. `timeout` is for what a long computation produces (a search's hits). */
export async function expectVisible(loc: Locator, visible: boolean, timeout?: number): Promise<void> {
  const shown = loc.filter({visible: true});
  await (visible ? expect(shown, 'visible expected').not.toHaveCount(0, {timeout}) : expect(shown, 'hidden expected').toHaveCount(0, {timeout}));
}

/** `aria-expanded` sits on the element itself (a tree row) or on its header/trigger inside; the
 * Dart tree says it with its twistie's class instead. `null` — the element says nothing at all —
 * is reported as that rather than as the opposite state. A u2 control's status says `expanded`, or
 * `open` for a popup. */
async function expectExpanded(on: On, loc: Locator, expanded: boolean): Promise<void> {
  if (await reports(loc, ['expanded', 'open']))
    return expectStatus(on, loc, expanded, ['expanded', 'collapsed'], (s) => (s.expanded ?? s.open) === true,
      (d) => d.expanded === null ? null : d.expanded === 'true', 'aria-expanded');
  await expect.poll(() => readExpanded(loc), {message: 'expanded state (aria-expanded, or the tree twistie)'})
    .toBe(expanded);
}

/** A branch whose children disagree says `aria-checked="mixed"` (a hierarchical filter node). */
async function expectMixed(on: On, loc: Locator, mixed: boolean): Promise<void> {
  if (await reports(loc, ['checked']))
    return expectStatus(on, loc, mixed, ['partially checked', 'not partially checked'], (s) => s.checked === 'mixed',
      (d) => d.checked === null ? null : d.checked === 'mixed', 'aria-checked');
  const box = loc.locator('[aria-checked]').first();
  const target = await box.count() > 0 ? box : loc.first();
  await (mixed ? expect(target) : expect(target).not).toHaveAttribute('aria-checked', 'mixed');
}

async function expectChecked(on: On, loc: Locator, checked: boolean): Promise<void> {
  await whenPresent(loc);
  if (await reports(loc, ['checked']))
    return expectStatus(on, loc, checked, ['checked', 'unchecked'], (s) => s.checked === true,
      (d) => d.checked === null ? null : d.checked === 'true', 'aria-checked');
  const box = loc.locator(CHECKABLE).first();
  const target = await box.count() > 0 ? box : loc;
  await (checked ? expect(target) : expect(target).not).toBeChecked();
}

/** Several matches (stacked notifications, repeated rows) mean "any of them" for a positive check
 * and "none of them" for a negative one, read over the visible ones when there are any: a hidden
 * copy (a tab panel not shown keeps its own) answers only where nothing is visible. A single match
 * is checked as itself. A tooltip is what is shown now: the platform keeps one tooltip element,
 * hidden between hovers. */
export async function expectText(page: Page, target: ElementRef, text: string, options: {exact?: boolean; negate?: boolean} = {}): Promise<void> {
  const plan = refOf(page, target).plan;
  const tooltip = plan.type === 'kind' && plan.kind.name === 'tooltip';
  const loc = await locate(page, target);
  if (tooltip || await loc.count() > 1) {
    const shown = loc.filter({visible: true});
    const want = options.exact ? exactText(text) : new RegExp(escapeRegExp(text), 'i');
    const hits = async () => (tooltip || await shown.count() > 0 ? shown : loc).filter({hasText: want}).count();
    const message = `the ${tooltip ? 'visible ' : ''}matches of ${target.phrase} with the text "${text}"`;
    await (options.negate ? expect.poll(hits, {message}).toBe(0) : expect.poll(hits, {message}).not.toBe(0));
    return;
  }
  const e = options.negate ? expect(loc).not : expect(loc);
  if (options.exact)
    await e.toHaveText(exactText(text));
  else
    await e.toContainText(text, {ignoreCase: true});
}

function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

/** A number the reader sees, between two bounds: what a slider or a stepper arrives at, where no
 * exact value is the claim. */
export async function expectValueBetween(page: Page, target: ElementRef, lo: number, hi: number): Promise<void> {
  await expect.poll(async () => {
    const text = await readValue(page, target);
    const n = text === undefined || text === '' ? NaN : Number(text);
    return n >= lo && n <= hi ? 'in range' : `${text ?? 'nothing'}`;
  }, {message: `the value of ${target.phrase}, expected between ${lo} and ${hi}`}).toBe('in range');
}

/** The value a reader would see. The editor is resolved on every attempt, not once: a view the
 * platform is still building exposes the input host before the control inside it, and a single
 * resolution then reads the `div.ui-input-editor` wrapper and fails with "Not an input element"
 * instead of waiting for the input that is about to arrive (a model reopened from its own address).
 * `undefined` is "nothing can answer yet", which keeps the poll going. A `<select>` answers with the
 * text of the chosen option, which is what the reader sees. A u2 control root answers with its status:
 * the text it shows, else its value. */
export async function readValue(page: Page, target: ElementRef): Promise<string | undefined> {
  return (await valueOf(page, target)).value;
}

async function valueOf(page: Page, target: ElementRef): Promise<{value: string | undefined; state?: U2State}> {
  const all = await locate(page, target).then((loc) => loc.evaluateAll(readControls, {errorPart: null, find: null})).catch(() => []);
  const r = shownFirst(all);
  const text = r?.state ? textOf(r.state) : undefined;
  if (r?.state && text !== undefined) {
    if (r.dom.value !== null && r.dom.value !== text)
      noteContract(page, target.phrase, `"${text}"`, `"${r.dom.value}"`, 'the editor\'s value');
    return {value: text, state: r.state};
  }
  // `editorOf` reads the first match's editor, so a hidden first match leaves the reading to the visible one
  if (r !== null && r !== all[0] && r.dom.value !== null)
    return {value: r.dom.value};
  const editor = await editorOf(page, target).catch(() => undefined);
  if (editor === undefined)
    return {value: undefined};
  return {value: await editor.first().evaluate((e) => {
    if (e instanceof HTMLSelectElement)
      return e.selectedIndex < 0 ? undefined : (e.options[e.selectedIndex].textContent ?? '').trim();
    if (e instanceof HTMLInputElement || e instanceof HTMLTextAreaElement)
      return e.value;
    return (e as HTMLElement).isContentEditable ? (e.textContent ?? '') : undefined;
  }).catch(() => undefined)};
}

/** Fails with `the value of Stage input is "Discovery", expected "Phase II" (u2: choice-input "Stage",
 * enabled, valid)`; when nothing could hold a value, as the wait it was. */
export async function expectValue(page: Page, target: ElementRef, value: string, negate = false): Promise<void> {
  const message = `the value of ${target.phrase}`;
  let last: {value: string | undefined; state?: U2State} = {value: undefined};
  const read = async () => (last = await valueOf(page, target)).value;
  try {
    if (!negate) {
      await expect.poll(read, {message}).toBe(value);
      return;
    }
    await expect.poll(read, {message: `${message} (something that can hold one)`}).not.toBe(undefined);
    await expect.poll(read, {message}).not.toBe(value);
  }
  catch (e) {
    if (last.value === undefined)
      throw e;
    throw new Error(`${message} is "${last.value}", expected ${negate ? 'anything but ' : ''}"${value}"` +
      (last.state ? ` (u2: ${summary(last.state, false)})` : ''));
  }
}

/** The choices a dropdown offers, in order: a native `<select>`'s options, else the element's own
 * option rows. The blank option of a nullable dropdown is no choice, and a phrase cannot name it. */
function dropdownChoices(loc: Locator): Promise<string[]> {
  return loc.first().evaluate((e) => {
    const select = e.tagName === 'SELECT' ? e as HTMLSelectElement : e.querySelector('select');
    const items = select ? Array.from(select.options) : Array.from(e.querySelectorAll('[role="option"]'));
    return items.map((o) => (o.textContent ?? '').trim()).filter(Boolean);
  }).catch(() => [] as string[]);
}

export async function expectOptions(page: Page, target: ElementRef, list: string): Promise<void> {
  const want = list.split(/\s*,\s*/).filter(Boolean);
  const loc = await locate(page, target);
  await openPropertyEditor(loc);
  await expect.poll(() => dropdownChoices(loc), {message: `the choices ${target.phrase} offers`}).toEqual(want);
}

/** One choice among those a dropdown offers, or missing from them; the dropdown must offer something,
 * so an empty or missing one cannot pass. */
export async function expectOption(page: Page, target: ElementRef, option: string, offered: boolean): Promise<void> {
  const loc = await locate(page, target);
  await openPropertyEditor(loc);
  const message = `the choices ${target.phrase} offers`;
  if (offered) {
    await expect.poll(() => dropdownChoices(loc), {message}).toContain(option);
    return;
  }
  await expect.poll(async () => (await dropdownChoices(loc)).length, {message: `${message}: any`}).toBeGreaterThan(0);
  await expect.poll(() => dropdownChoices(loc), {message}).not.toContain(option);
}

/** Rows of a collection: the first row vocabulary that has any is the one counted. */
export async function expectCount(page: Page, target: ElementRef, count: number): Promise<void> {
  const loc = await locate(page, target);
  for (const rows of ROWS) {
    const items = loc.locator(rows);
    if (await items.count() > 0)
      return expect(items).toHaveCount(count);
  }
  await expect(loc.locator(ROWS.join(', '))).toHaveCount(count);
}
