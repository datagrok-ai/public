/* What a u2 control says about itself (`DG.U2.state`, the u2 ↔ bdd contract §1.1), read for the generic
   steps — only when the located element is the control's root: a part, a wrapper and a Dart widget keep
   their DOM readings, and a key the control does not report falls back to the DOM too (absence is not a
   value). When the status and the DOM both answer and disagree, the claim follows the status and the
   disagreement is printed once per phrase and page, a finding for the u2 side. */
import type {Locator, Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {describeCandidates, locate, refOf} from './locate.js';
import {expect} from './patience.js';
import type {NounRef} from '../nouns.js';

export type JsonValue = null | boolean | number | string | JsonValue[] | {[key: string]: JsonValue};

export interface U2ItemState {
  key: string;
  label: string;
  index: number;
  selected?: boolean;
  checked?: boolean | 'mixed';
  expanded?: boolean;
  level?: number;
  disabled?: boolean;
  current?: boolean;
  completed?: boolean;
  invalid?: boolean;
  rendered: boolean;
}

export interface U2ItemsState {
  count: number;
  list: U2ItemState[];
  truncated: boolean;
}

/** `IControlState` of `u2core/protocol.ts`: keys a kind does not have are absent. */
export interface U2State {
  kind: string;
  tag?: string;
  name?: string;
  label?: string;
  enabled: boolean;
  visible: boolean;
  focused: boolean;
  readOnly?: boolean;
  required?: boolean;
  valid?: boolean;
  error?: string | null;
  busy?: boolean;
  value?: JsonValue;
  text?: string;
  placeholder?: string;
  open?: boolean;
  expanded?: boolean;
  checked?: boolean | 'mixed';
  pressed?: boolean;
  selected?: boolean;
  items?: U2ItemsState;
  selection?: string[];
  current?: string;
  extra?: {[key: string]: number | string | boolean};
}

/** What the DOM says on the element or its focus target (`<rootId>-editor`, else the first editable
 * element inside): the ARIA states, and the text the editor shows (`null` where none is readable). */
export interface DomReadings {
  disabled: boolean;
  invalid: boolean;
  readOnly: boolean;
  busy: boolean;
  required: boolean;
  checked: string | null;
  expanded: string | null;
  value: string | null;
}

export interface ControlReading {
  shown: boolean;
  /** The status, when the element is a u2 control root. */
  state: U2State | null;
  dom: DomReadings;
  /** The text of the element's error part: `null` without one, '' while it is hidden or empty. */
  errorText: string | null;
  /** Items looked up by key, else by label (`find`), and the selected and current items (`picked`). */
  found: {matches: U2ItemState[][]; picked: U2ItemState[]} | null;
}

export interface ReadOptions {
  /** The selector of the error part. */
  errorPart?: string | null;
  /** Item keys or labels to look up among all the items. */
  find?: string[] | null;
}

/** In the page: one reading per element. */
export const readControls = (els: Element[], o: {errorPart: string | null; find: string[] | null}): ControlReading[] => {
  const U2 = (window as any).DG?.U2;
  return els.map((el) => {
    let state: any = null;
    let found: ControlReading['found'] = null;
    try {
      state = U2?.control?.(el)?.root === el ? U2.state(el) : null;
      if (state && o.find) {
        const all: any[] = U2.items?.(el, {all: true})?.list ?? state.items?.list ?? [];
        const keyed = new Set<string>([...(state.selection ?? []), ...(state.current != null ? [state.current] : [])]);
        found = {
          matches: o.find.map((w) => {
            const byKey = all.filter((i) => i.key === w);
            return byKey.length > 0 ? byKey : all.filter((i) => String(i.label).trim().toLowerCase() === w.trim().toLowerCase());
          }),
          picked: all.filter((i) => keyed.has(i.key)),
        };
      }
    }
    catch {
      state = null;
      found = null;
    }
    const target = (el.id ? document.getElementById(`${el.id}-editor`) : null) ?? el.querySelector('input, select, textarea, ' +
      '[role="combobox"], [role="listbox"], [role="slider"], [role="switch"], [role="checkbox"], [contenteditable="true"]') ?? el;
    const aria: Record<string, string | null> = Object.fromEntries(['aria-disabled', 'aria-invalid', 'aria-readonly', 'aria-required',
      'aria-checked', 'aria-expanded'].map((n) => [n, el.getAttribute(n) ?? (target !== el ? target.getAttribute(n) : null)]));
    const input = target as HTMLInputElement;
    const box = target instanceof HTMLInputElement && (target.type === 'checkbox' || target.type === 'radio');
    const value = box ? null : target instanceof HTMLSelectElement ?
      (target.selectedIndex < 0 ? '' : (target.options[target.selectedIndex].textContent ?? '').trim()) :
      target instanceof HTMLInputElement || target instanceof HTMLTextAreaElement ? target.value :
        (target as HTMLElement).isContentEditable ? (target.textContent ?? '') : null;
    const error = o.errorPart ? el.querySelector(o.errorPart) : null;
    return {
      shown: el.getClientRects().length > 0 && getComputedStyle(el).visibility !== 'hidden',
      state,
      dom: {
        disabled: aria['aria-disabled'] === 'true' || input.disabled === true,
        invalid: aria['aria-invalid'] === 'true',
        readOnly: aria['aria-readonly'] === 'true' || input.readOnly === true,
        busy: el.closest('[aria-busy="true"]') !== null,
        required: aria['aria-required'] === 'true' || input.required === true,
        checked: box ? (input.indeterminate ? 'mixed' : String(input.checked)) : aria['aria-checked'],
        expanded: aria['aria-expanded'],
        value,
      },
      errorText: error === null ? null : error.getClientRects().length === 0 ? '' : (error.textContent ?? '').replace(/\s+/g, ' ').trim(),
      found,
    };
  });
};

/** Of readings in DOM order, the first shown one, else the first: a hidden copy (a tab panel not shown
 * keeps its own) answers only where nothing is shown. */
export function shownFirst<T extends {shown: boolean}>(all: T[]): T | null {
  return all.find((r) => r.shown) ?? all[0] ?? null;
}

/** The first visible element of `loc`, else its first, read in one roundtrip without waiting for it;
 * `null` when there is none. */
export async function readState(loc: Locator, options: ReadOptions = {}): Promise<ControlReading | null> {
  return shownFirst(await loc.evaluateAll(readControls, {errorPart: options.errorPart ?? null, find: options.find ?? null}));
}

/** The items of the u2 control whose root is the first element of `loc` (uncapped with `all`); `null`
 * for anything else. */
export async function readItems(loc: Locator, all = false): Promise<U2ItemsState | null> {
  return (await loc.first().evaluateAll((els, everything) => els.map((el) => {
    const U2 = (window as any).DG?.U2;
    return U2?.control?.(el)?.root === el ? U2.items(el, {all: everything}) ?? null : null;
  }), all))[0] ?? null;
}

/** Whether the status of the first element of `loc` (with `single`: of its only element) has one of `keys`. */
export async function reports(loc: Locator, keys: (keyof U2State)[], single = false): Promise<boolean> {
  const all = await (single ? loc : loc.first()).evaluateAll(readControls, {errorPart: null, find: null});
  const state = all.length === 1 ? all[0].state : null;
  return state !== null && keys.some((k) => k in state);
}

/** The value as the user reads it: the displayed text, else the raw value (an array comma-joined). */
export function textOf(s: U2State): string | undefined {
  if (typeof s.text === 'string')
    return s.text;
  if (s.value === undefined)
    return undefined;
  return Array.isArray(s.value) ? s.value.map((v) => typeof v === 'object' && v !== null ? JSON.stringify(v) : String(v)).join(', ') :
    s.value === null ? '' : typeof s.value === 'object' ? JSON.stringify(s.value) : String(s.value);
}

/** One line for a failure: `choice-input "Stage", enabled, valid, value "Discovery"`. */
export function summary(s: U2State, withValue = true): string {
  const parts = [`${s.kind}${s.label ? ` "${s.label}"` : ''}`, s.enabled ? 'enabled' : 'disabled'];
  if (s.visible === false)
    parts.push('hidden');
  for (const [key, word] of [['readOnly', 'read-only'], ['required', 'required'], ['busy', 'busy']] as const) {
    if (s[key] === true)
      parts.push(word);
  }
  if (s.valid !== undefined)
    parts.push(s.valid ? 'valid' : s.error ? `invalid ("${s.error}")` : 'invalid');
  const text = withValue ? textOf(s) : undefined;
  if (text !== undefined)
    parts.push(`value "${text}"`);
  return parts.join(', ');
}

const noted = new WeakMap<Page, Set<string>>();

/** A status/DOM disagreement, printed once per phrase, attribute and page. */
export function noteContract(page: Page, phrase: string, status: string, dom: string, attribute: string): void {
  if (!noted.has(page))
    noted.set(page, new Set());
  const seen = noted.get(page)!;
  if (seen.has(`${phrase}\n${attribute}`))
    return;
  seen.add(`${phrase}\n${attribute}`);
  console.warn(`bdd: u2 contract: ${phrase}: the status says ${status}, the DOM says ${dom} (${attribute})`);
}

/** The sentence a status-only claim fails with on an element that has no status, or no such key. */
export function reportsNo(phrase: string, r: ControlReading, what: string): string {
  return r.state === null ? `${phrase} is not a u2 control; it reports no ${what}` : `${phrase} (u2: ${summary(r.state)}) reports no ${what}`;
}

function nounOf(ref: NounRef): string {
  return ref.plan.type === 'kind' ? ref.plan.kind.name : ref.plan.type === 'part' ? ref.plan.part : 'element';
}

/** A claim of a u2 phrase, polled: the visible matches are read in one roundtrip each time, and `judge`
 * says '' when the claim holds, else why not in one sentence (a judge that throws ends the poll with
 * its message). Nothing visible fails as a wait, so the report lists what the page shows; several
 * visible matches without an ordinal fail naming them (R9). */
export async function expectControl(page: Page, target: ElementRef, judge: (r: ControlReading) => string,
  options: ReadOptions = {}): Promise<void> {
  let reason = '';
  let ambiguous: Locator | undefined;
  try {
    await expect.poll(async () => {
      const ref = refOf(page, target);
      const loc = (await locate(page, target)).filter({visible: true});
      const all = await loc.evaluateAll(readControls, {errorPart: options.errorPart ?? null, find: options.find ?? null});
      ambiguous = all.length > 1 && ref.ordinal === undefined ? loc : undefined;
      try {
        reason = all.length === 0 ? '' : ambiguous ? `"${target.phrase}" matches ${all.length} visible ${nounOf(ref)}s` : judge(all[0]);
      }
      catch (e) {
        reason = e instanceof Error ? e.message : String(e);
        throw e;
      }
      return all.length > 0 && reason === '';
    }, {message: `${target.phrase} as a u2 control`}).toBe(true);
  }
  catch (e) {
    if (!reason)
      throw e;
    throw new Error(ambiguous ? `${reason}: ${await describeCandidates(ambiguous).catch(() => 'unreadable')}` : reason);
  }
}
