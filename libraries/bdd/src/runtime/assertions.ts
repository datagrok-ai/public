/* Outcome checks over Playwright's retrying `expect`, so a feature never needs a wait step. */
import {expect, Locator, Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {editorOf} from './gestures.js';
import {exactText, locate} from './locate.js';

export type State = 'visible' | 'hidden' | 'present' | 'absent' | 'enabled' | 'disabled' | 'checked' |
  'unchecked' | 'selected' | 'empty' | 'expanded' | 'collapsed' | 'focused';

export const STATES: State[] = ['visible', 'hidden', 'present', 'absent', 'enabled', 'disabled', 'checked',
  'unchecked', 'selected', 'empty', 'expanded', 'collapsed', 'focused'];

const ROWS = ['.u2-list-row', '[role="option"]', '[role="row"]', '[role="tab"]', 'option', '.d4-list-item', 'tbody tr', 'tr', 'li'];
const SELECTED = '[aria-selected="true"], [aria-pressed="true"], [aria-checked="true"], [aria-current]:not([aria-current="false"]), ' +
  '.u2-list-row-selected';

async function control(page: Page, target: ElementRef, loc: Locator): Promise<Locator> {
  const inner = loc.locator('input, select, textarea, button').first();
  return await inner.count() > 0 ? await editorOf(page, target) : loc;
}

export async function expectState(page: Page, target: ElementRef, state: State, negate = false): Promise<void> {
  const loc = await locate(page, target);
  const e = negate ? expect(loc).not : expect(loc);
  switch (state) {
    case 'visible': return expectVisible(loc, !negate);
    case 'hidden': return expectVisible(loc, negate);
    case 'present': return negate ? expect(loc).toHaveCount(0) : expect(loc.first()).toBeAttached();
    case 'absent': return negate ? expect(loc.first()).toBeAttached() : expect(loc).toHaveCount(0);
    case 'enabled': return expectEnabled(await control(page, target, loc), !negate);
    case 'disabled': return expectEnabled(await control(page, target, loc), negate);
    case 'checked': return expectChecked(loc, !negate);
    case 'unchecked': return expectChecked(loc, negate);
    case 'selected': return expectSelected(page, loc, !negate);
    case 'empty': return (negate ? expect(await editorOf(page, target)).not : expect(await editorOf(page, target))).toHaveValue('');
    case 'expanded': return expectExpanded(loc, !negate);
    case 'collapsed': return expectExpanded(loc, negate);
    case 'focused': return e.toBeFocused();
  }
}

/** Several matches (stacked balloons, repeated rows): visible when any is, hidden when none is. */
async function expectVisible(loc: Locator, visible: boolean): Promise<void> {
  if (await loc.count() > 1) {
    const shown = loc.filter({visible: true});
    await (visible ? expect(shown).not.toHaveCount(0) : expect(shown).toHaveCount(0));
    return;
  }
  await (visible ? expect(loc).toBeVisible() : expect(loc).toBeHidden());
}

/** Options and tabs say `aria-selected`, toggles and cards `aria-pressed`, radio-like buttons
 * `aria-checked`, wizard steps and breadcrumbs `aria-current`. */
async function expectSelected(page: Page, loc: Locator, selected: boolean): Promise<void> {
  const hit = loc.and(page.locator(SELECTED));
  await (selected ? expect(hit).not.toHaveCount(0) : expect(hit).toHaveCount(0));
}

/** `aria-expanded` sits on the element itself (a tree row) or on its header/trigger inside. */
async function expectExpanded(loc: Locator, expanded: boolean): Promise<void> {
  const self = loc.first();
  const own = await self.getAttribute('aria-expanded', {timeout: 2000}).catch(() => null);
  const control = own !== null ? self : self.locator('[aria-expanded]').first();
  await expect(control).toHaveAttribute('aria-expanded', String(expanded));
}

async function expectEnabled(loc: Locator, enabled: boolean): Promise<void> {
  if (enabled) {
    await expect(loc).toBeEnabled();
    await expect(loc).not.toHaveClass(/u2-input-disabled|d4-disabled/);
  }
  else
    await expect(loc.or(loc.locator('xpath=ancestor-or-self::*[contains(@class, "u2-input-disabled")]'))).toBeDisabled().catch(() =>
      expect(loc).toHaveClass(/u2-input-disabled|d4-disabled/));
}

async function expectChecked(loc: Locator, checked: boolean): Promise<void> {
  const box = loc.locator('input[type="checkbox"], input[type="radio"], [role="checkbox"], [role="switch"]').first();
  const target = await box.count() > 0 ? box : loc;
  await (checked ? expect(target) : expect(target).not).toBeChecked();
}

/** Several matches (stacked notifications, repeated rows) mean "any of them" for a positive check
 * and "none of them" for a negative one; a single match is checked as itself. */
export async function expectText(page: Page, target: ElementRef, text: string, options: {exact?: boolean; negate?: boolean} = {}): Promise<void> {
  const loc = await locate(page, target);
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
