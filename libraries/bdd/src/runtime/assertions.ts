/* Outcome checks over Playwright's retrying `expect`, so a feature never needs a wait step. */
import {Locator, Page} from '@playwright/test';
import {expect} from './patience.js';
import type {ElementRef} from './args.js';
import {editorOf, readExpanded} from './gestures.js';
import {exactText, locate, locateActionable, refOf} from './locate.js';

import type {State} from '../states.js';

export type {State};
export {STATES} from '../states.js';
const INVALID_CLASSES = ['d4-invalid', 'd4-forced-invalid', 'u2-input-invalid'];

const ROWS = ['.u2-list-row', '[role="option"]', '[role="row"]', '[role="tab"]', 'option', '.d4-list-item', '[name="legend-item"]', 'tbody tr', 'tr', 'li'];
const SELECTED = '[aria-selected="true"], [aria-pressed="true"], [aria-checked="true"], [aria-current]:not([aria-current="false"]), ' +
  '.u2-list-row-selected';

export async function expectState(page: Page, target: ElementRef, state: State, negate = false): Promise<void> {
  const loc = ['visible', 'hidden', 'present', 'absent', 'enabled', 'disabled'].includes(state) ?
    await locate(page, target) : await locateActionable(page, target);
  const e = negate ? expect(loc).not : expect(loc);
  switch (state) {
    case 'visible': return expectVisible(loc, !negate);
    case 'hidden': return expectVisible(loc, negate);
    case 'present': return negate ? expect(loc).toHaveCount(0) : expect(loc.first()).toBeAttached();
    case 'absent': return negate ? expect(loc.first()).toBeAttached() : expect(loc).toHaveCount(0);
    case 'enabled': return expectEnabled(loc, !negate);
    case 'disabled': return expectEnabled(loc, negate);
    case 'checked': return expectChecked(loc, !negate);
    case 'unchecked': return expectChecked(loc, negate);
    case 'partially checked': return expectMixed(loc, !negate);
    case 'invalid': return expectInvalid(loc, !negate);
    case 'valid': return expectInvalid(loc, negate);
    case 'selected': return expectSelected(page, loc, !negate);
    case 'empty': return (negate ? expect(await editorOf(page, target)).not : expect(await editorOf(page, target))).toHaveValue('');
    case 'expanded': return expectExpanded(loc, !negate);
    case 'collapsed': return expectExpanded(loc, negate);
    case 'focused': return e.toBeFocused();
  }
}

/** Several matches (stacked balloons, repeated rows): visible when any is, hidden when none is —
 * one query either way. */
async function expectVisible(loc: Locator, visible: boolean): Promise<void> {
  const shown = loc.filter({visible: true});
  await (visible ? expect(shown, 'visible expected').not.toHaveCount(0) : expect(shown, 'hidden expected').toHaveCount(0));
}

/** Options and tabs say `aria-selected`, toggles and cards `aria-pressed`, radio-like buttons
 * `aria-checked`, wizard steps and breadcrumbs `aria-current`. */
async function expectSelected(page: Page, loc: Locator, selected: boolean): Promise<void> {
  const hit = loc.and(page.locator(SELECTED));
  await (selected ? expect(hit).not.toHaveCount(0) : expect(hit).toHaveCount(0));
}

/** `aria-expanded` sits on the element itself (a tree row) or on its header/trigger inside; the
 * Dart tree says it with its twistie's class instead. `null` — the element says nothing at all —
 * is reported as that rather than as the opposite state. */
async function expectExpanded(loc: Locator, expanded: boolean): Promise<void> {
  await expect.poll(() => readExpanded(loc), {message: 'expanded state (aria-expanded, or the tree twistie)'})
    .toBe(expanded);
}

/** Disabled: the element or an ancestor says so (`aria-disabled`, the u2/Dart disabled classes),
 * or it — or the control inside it — is natively disabled. Over the visible matches when there
 * are any (the Dart menu's hidden mirror), else all of them (a property row in a panel that is
 * not shown still says whether it is gated). One query per poll. */
async function expectEnabled(loc: Locator, enabled: boolean): Promise<void> {
  const disabled = () => loc.evaluateAll((all) => {
    const shown = all.filter((e) => e.getClientRects().length > 0 && getComputedStyle(e).visibility !== 'hidden');
    const els = shown.length > 0 ? shown : all;
    if (els.length === 0)
      return undefined;
    const marked = (e: Element) => e.getAttribute('aria-disabled') === 'true' ||
      ['u2-input-disabled', 'd4-disabled', 'd4-menu-item-disabled'].some((c) => e.classList.contains(c));
    return els.every((el) => {
      for (let e: Element | null = el; e; e = e.parentElement) {
        if (marked(e))
          return true;
      }
      const native = el.matches('input, select, textarea, button, fieldset, option') ? el : el.querySelector('input, select, textarea, button');
      return native !== null && (native as HTMLInputElement).disabled === true;
    });
  });
  await expect.poll(disabled, {message: `${enabled ? 'enabled' : 'disabled'} expected`}).toBe(!enabled);
}

/** A branch whose children disagree says `aria-checked="mixed"` (a hierarchical filter node). */
async function expectMixed(loc: Locator, mixed: boolean): Promise<void> {
  const box = loc.locator('[aria-checked]').first();
  const target = await box.count() > 0 ? box : loc.first();
  await (mixed ? expect(target) : expect(target).not).toHaveAttribute('aria-checked', 'mixed');
}

/** Invalid: `aria-invalid` or the platform's invalid classes on the element or a control inside. */
async function expectInvalid(loc: Locator, invalid: boolean): Promise<void> {
  const holds = () => loc.evaluateAll((all, classes: string[]) => all.some((el) => {
    const marked = (e: Element) => e.getAttribute('aria-invalid') === 'true' || classes.some((c) => e.classList.contains(c));
    return marked(el) || Array.from(el.querySelectorAll('input, select, textarea, [aria-invalid]')).some(marked);
  }), INVALID_CLASSES);
  await expect.poll(holds, {message: `${invalid ? 'invalid' : 'valid'} expected`}).toBe(invalid);
}

async function expectChecked(loc: Locator, checked: boolean): Promise<void> {
  const box = loc.locator('input[type="checkbox"], input[type="radio"], [role="checkbox"], [role="switch"]').first();
  const target = await box.count() > 0 ? box : loc;
  await (checked ? expect(target) : expect(target).not).toBeChecked();
}

/** Several matches (stacked notifications, repeated rows) mean "any of them" for a positive check
 * and "none of them" for a negative one; a single match is checked as itself. A tooltip is what
 * is shown now: the platform keeps one tooltip element, hidden between hovers. */
export async function expectText(page: Page, target: ElementRef, text: string, options: {exact?: boolean; negate?: boolean} = {}): Promise<void> {
  const plan = refOf(page, target).plan;
  const loc = plan.type === 'kind' && plan.kind.name === 'tooltip' ? await locateActionable(page, target) : await locate(page, target);
  if (await loc.count() > 1) {
    const matching = loc.filter({hasText: options.exact ? exactText(text) : new RegExp(escapeRegExp(text), 'i')});
    await (options.negate ? expect(matching).toHaveCount(0) : expect(matching).not.toHaveCount(0));
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

export async function expectValue(page: Page, target: ElementRef, value: string): Promise<void> {
  const editor = await editorOf(page, target);
  const tag = await editor.evaluate((e) => e.tagName).catch(() => '');
  if (tag === 'SELECT') {
    const selected = editor.locator('option:checked');
    await expect(selected).toHaveText(exactText(value));
  }
  else
    await expect(editor).toHaveValue(value);
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
