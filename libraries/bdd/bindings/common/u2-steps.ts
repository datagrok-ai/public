/* Claims a u2 control answers through its status (`DG.U2.state`, R10): its error, all of its items, an
   item's state, its selection, the current tab or step, its raw value. Steps only — no nouns. An element
   that is not a u2 control root fails saying so, unless the step names a DOM reading to fall back to;
   several visible matches of the phrase fail naming them (R9). */
import {Page} from '@playwright/test';
import {Then} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {refOf} from '../../src/runtime/locate.js';
import {ControlReading, expectControl, reportsNo, U2ItemState} from '../../src/runtime/u2-state.js';
import type {State} from '../../src/states.js';

const LISTED_MAX = 20;

/** The error the status reports, else the error part's text; `null` when neither exists. */
function errorOf(r: ControlReading): string | null {
  return r.state && 'error' in r.state ? r.state.error ?? '' : r.errorText;
}

function errorPart(page: Page, target: ElementRef): string {
  const plan = refOf(page, target).plan;
  return (plan.type === 'kind' ? plan.kind.parts?.error : undefined) ?? '[data-u2-part="error"]';
}

function listed(r: ControlReading): string {
  const items = r.state?.items;
  if (!items || items.count === 0)
    return 'nothing';
  return items.list.slice(0, LISTED_MAX).map((i) => `"${i.label}"`).join(', ') + (items.count > LISTED_MAX ? `, … ${items.count - LISTED_MAX} more` : '');
}

function selectionOf(r: ControlReading): string {
  const keys = r.state?.selection ?? [];
  return keys.length === 0 ? 'nothing' : `"${keys.map((k) => r.found?.picked.find((i) => i.key === k)?.label ?? k).join(', ')}"`;
}

function fieldsOf(item: U2ItemState): string {
  return Object.entries(item).filter(([k]) => !['key', 'label', 'index'].includes(k)).map(([k, v]) => `${k} ${v}`).join(', ');
}

export const shouldHaveError = Then('{element} should have the error {string}', (page: Page, target: ElementRef, error: string) =>
  expectControl(page, target, (r) => {
    const shown = errorOf(r);
    if (shown === null)
      return reportsNo(target.phrase, r, 'error');
    return shown === error ? '' : `the error of ${target.phrase} is ${shown ? `"${shown}"` : 'none'}, expected "${error}"`;
  }, {errorPart: errorPart(page, target)}),
{description: 'the validation message the u2 status reports, exactly; else the text of the error part'});

export const shouldHaveNoError = Then('{element} should have no error', (page: Page, target: ElementRef) =>
  expectControl(page, target, (r) => {
    const shown = errorOf(r);
    if (shown === null)
      return reportsNo(target.phrase, r, 'error');
    return shown === '' ? '' : `${target.phrase} has the error "${shown}"`;
  }, {errorPart: errorPart(page, target)}),
{description: 'the u2 status reports no validation message; else the error part is hidden or empty'});

export const shouldHaveItemsInTotal = Then('{element} should have {int} item(s) in total', (page: Page, target: ElementRef, count: number) =>
  expectControl(page, target, (r) => {
    const items = r.state?.items;
    if (!items)
      return reportsNo(target.phrase, r, 'items');
    return items.count === count ? '' : `${target.phrase} has ${items.count} item${items.count === 1 ? '' : 's'} in total, expected ${count}`;
  }),
{description: 'every item the u2 control holds, the rows a virtual list has not rendered included; "should have N items" counts the rendered rows'});

export const shouldListItem = Then('{element} should list the item {string}', (page: Page, target: ElementRef, item: string) =>
  expectControl(page, target, (r) => {
    if (!r.state?.items || !r.found)
      return reportsNo(target.phrase, r, 'items');
    return r.found.matches[0].length > 0 ? '' : `${target.phrase} does not list "${item}"; it lists: ${listed(r)}`;
  }, {find: [item]}),
{description: 'an item of the u2 control, rendered or not, by its key, else by its label (any case)'});

export const shouldNotListItem = Then('{element} should not list the item {string}', (page: Page, target: ElementRef, item: string) =>
  expectControl(page, target, (r) => {
    if (!r.state?.items || !r.found)
      return reportsNo(target.phrase, r, 'items');
    const [match] = r.found.matches;
    return match.length === 0 ? '' : `${target.phrase} lists "${item}" (${match.map((i) => `key ${i.key}`).join(', ')})`;
  }, {find: [item]}),
{description: 'no item of the u2 control, rendered or not, has that key or label (any case)'});

/** The key an item reports a state under and the reading; disabled and invalid are reported only when
 * true, as their ARIA is, so their absence reads as enabled and valid. */
const ITEM_STATES: {[s in State]?: [keyof U2ItemState, (v: unknown) => boolean]} = {
  'selected': ['selected', (v) => v === true],
  'checked': ['checked', (v) => v === true],
  'unchecked': ['checked', (v) => v === false],
  'partially checked': ['checked', (v) => v === 'mixed'],
  'expanded': ['expanded', (v) => v === true],
  'collapsed': ['expanded', (v) => v === false],
  'disabled': ['disabled', (v) => v === true],
  'enabled': ['disabled', (v) => v !== true],
  'invalid': ['invalid', (v) => v === true],
  'valid': ['invalid', (v) => v !== true],
};

export const itemShouldBe = Then('item {string} of {element} should be {state}', (page: Page, item: string, target: ElementRef, state: State) =>
  expectControl(page, target, (r) => {
    if (!r.state?.items || !r.found)
      return reportsNo(target.phrase, r, 'items');
    const [match] = r.found.matches;
    if (match.length === 0)
      return `${target.phrase} has no item "${item}"; it lists: ${listed(r)}`;
    if (match.length > 1)
      return `"${item}" matches ${match.length} items of ${target.phrase}: ${match.map((i) => i.key).join(', ')}`;
    const it = match[0];
    const reading = ITEM_STATES[state];
    if (!reading || (!(reading[0] in it) && reading[0] !== 'disabled' && reading[0] !== 'invalid'))
      throw new Error(`item "${item}" of ${target.phrase} reports no ${state} state; it reports: ${fieldsOf(it) || 'nothing'}`);
    return reading[1](it[reading[0]]) ? '' : `item "${item}" of ${target.phrase} is not ${state}; it reports: ${fieldsOf(it)}`;
  }, {find: [item]}),
{description: 'an item of the u2 control by key, else label: selected, (un/partially) checked, expanded or collapsed, enabled or disabled, valid or invalid'});

export const shouldHaveSelection = Then('{element} should have the selection {string}', (page: Page, target: ElementRef, list: string) => {
  const want = list.split(/\s*,\s*/).filter(Boolean);
  return expectControl(page, target, (r) => {
    if (!r.state?.selection || !r.found)
      return reportsNo(target.phrase, r, 'selection');
    const wanted = new Set<string>();
    for (const [i, match] of r.found.matches.entries()) {
      if (match.length > 1)
        return `"${want[i]}" matches ${match.length} items of ${target.phrase}: ${match.map((m) => m.key).join(', ')}`;
      wanted.add(match.length === 1 ? match[0].key : want[i]);
    }
    const keys = new Set(r.state.selection);
    return wanted.size === keys.size && [...wanted].every((k) => keys.has(k)) ? '' :
      `the selection of ${target.phrase} is ${selectionOf(r)}, expected "${want.join(', ')}"`;
  }, {find: want});
}, {description: 'the selected items of the u2 control as a set: a comma-separated list of keys or labels'});

export const shouldHaveNothingSelected = Then('{element} should have nothing selected', (page: Page, target: ElementRef) =>
  expectControl(page, target, (r) => {
    if (!r.state?.selection)
      return reportsNo(target.phrase, r, 'selection');
    return r.state.selection.length === 0 ? '' : `${target.phrase} has the selection ${selectionOf(r)}`;
  }, {find: []}));

function expectCurrent(page: Page, target: ElementRef, want: string, what: 'step' | 'tab'): Promise<void> {
  return expectControl(page, target, (r) => {
    if (!r.state || r.state.current === undefined || !r.found)
      return reportsNo(target.phrase, r, `current ${what}`);
    const key = r.state.current;
    const label = r.found.picked.find((i) => i.key === key)?.label ?? key;
    if (key === want || label.trim().toLowerCase() === want.trim().toLowerCase())
      return '';
    return what === 'step' ? `${target.phrase} is on the "${label}" step, expected "${want}"` : `${target.phrase} shows the "${label}" tab, expected "${want}"`;
  }, {find: []});
}

export const shouldBeOnStep = Then('{element} should be on the {string} step', (page: Page, target: ElementRef, step: string) =>
  expectCurrent(page, target, step, 'step'), {description: 'the current step of the u2 control by key, else by title (any case)'});

export const shouldShowTab = Then('{element} should show the {string} tab', (page: Page, target: ElementRef, tab: string) =>
  expectCurrent(page, target, tab, 'tab'), {description: 'the current tab of the u2 control by key, else by caption (any case)'});

export const shouldReportValue = Then('{element} should report the value {string}', (page: Page, target: ElementRef, value: string) => {
  let parsed: unknown;
  try {
    parsed = JSON.parse(value);
  }
  catch {
    parsed = undefined;
  }
  return expectControl(page, target, (r) => {
    if (!r.state || !('value' in r.state))
      return reportsNo(target.phrase, r, 'value');
    const v = r.state.value;
    return (parsed !== undefined && JSON.stringify(v) === JSON.stringify(parsed)) || String(v) === value ? '' :
      `${target.phrase} reports the value ${JSON.stringify(v)}, expected ${value}`;
  });
}, {description: 'the raw value of the u2 status, where the shown text formats it: compared as JSON when the argument parses as JSON, else as text'});
