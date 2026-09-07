import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';
import {cardCaptions} from '../../helpers/filter-panel';

const SEARCH = 'input.d4-column-selector-search-input';
const BACKDROP = '.d4-column-selector-backdrop';

export async function pickerSearchValue(page: Page): Promise<string> {
  return page.evaluate((s) => (document.querySelector(s) as HTMLInputElement | null)?.value ?? '', SEARCH);
}

// The picker builds its search box from the first keystroke and the value read straight after
// keyboard.type can still lag the last key ("RAC" for RACE), so the value is polled up to the
// typed text rather than read once.
export async function typeIntoOpenPicker(page: Page, column: string): Promise<string> {
  await page.keyboard.press(column[0].toLowerCase());
  await page.locator(SEARCH).first().waitFor({state: 'attached', timeout: 10_000});
  await page.keyboard.press('Control+a');
  await page.keyboard.type(column, {delay: 15});
  return v.pollValue(() => pickerSearchValue(page), (val) => val === column, 800, 30);
}

// Escape is only sent while the popup is up: with the panel focused it toggles the filter group.
export async function closePicker(page: Page): Promise<void> {
  if (await page.locator(BACKDROP).count() === 0) return;
  await page.keyboard.press('Escape');
  await expect.poll(() => page.locator(BACKDROP).count(), {
    timeout: 10_000,
    intervals: [30, 60, 120, 250, 500, 1000],
    message: 'the column picker popup stayed open and would intercept the next gesture',
  }).toBe(0);
}

// Enter in the open popup commits the row under the pointer when there is one and the typed
// name only otherwise (column_combo_box.dart:321 vs :361), so the pointer is parked away from
// the panel and the combo is opened with a dispatched mousedown rather than a real click.
export async function openHeaderPicker(page: Page): Promise<boolean> {
  await page.mouse.move(0, 0);
  const opened = await page.evaluate(() => {
    const combos = [...document.querySelectorAll(
      '[name="viewer-Filters"] .d4-filter-group-header [name="div-column-combobox-"]')];
    const combo = combos.find((e) => { const r = e.getBoundingClientRect(); return r.width > 0 && r.height > 0; });
    if (!combo) return false;
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const label = combo.querySelector('.d4-column-selector-column');
    (label ?? combo).dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
    return true;
  });
  expect(opened, 'the panel header exposes no laid-out column selector to add a filter with').toBe(true);
  return page.waitForFunction((b) => !!document.querySelector(b), BACKDROP, {timeout: 8000})
    .then(() => true).catch(() => false);
}

// Opens the laid-out panel's header picker and types the column until the search box provably
// holds it; a miss re-opens the picker once instead of committing a stray pick.
export async function typeColumnIntoHeaderPicker(page: Page, column: string): Promise<string> {
  let held = '';
  for (let attempt = 0; attempt < 2 && held !== column; attempt++) {
    if (attempt > 0) await closePicker(page);
    expect(await openHeaderPicker(page), 'the panel header column picker did not open').toBe(true);
    held = await typeIntoOpenPicker(page, column);
  }
  return held;
}

export async function addCardViaPicker(page: Page, column: string): Promise<void> {
  expect(await cardCaptions(page),
    `the panel must not already carry a ${column} card before the column selector adds it`)
    .not.toContain(column);
  let added = false;
  for (let attempt = 0; attempt < 2 && !added; attempt++) {
    expect(await typeColumnIntoHeaderPicker(page, column),
      `the column picker search box does not hold ${column}, so Enter would commit some other pick`)
      .toBe(column);
    await page.keyboard.press('Enter');
    added = await v.pollValue(async () => (await cardCaptions(page)).includes(column),
      (there) => there, attempt === 0 ? 5000 : 15_000, 50);
    if (!added) await closePicker(page);
  }
  expect(added, `no ${column} card came out of the column-picker commit`).toBe(true);
  // The caption lands before the card has finished growing, and a caller that measures a category
  // row off the fresh card would address a band the layout is about to move.
  await page.evaluate((col) => (window as any).__settledFor(() => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name') as HTMLElement | null)?.textContent ?? '').trim() === col);
    const r = card?.getBoundingClientRect();
    return r ? `${Math.round(r.top)}|${Math.round(r.height)}` : '';
  }, 150, 3000, 25), column);
  await expect.poll(() => page.locator(BACKDROP).count(), {
    timeout: 10_000,
    intervals: [30, 60, 120, 250, 500, 1000],
    message: 'the column picker popup stayed open and would intercept the next gesture',
  }).toBe(0);
}
