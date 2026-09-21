/* ---
realizes: [filters.cp.panel-core-ladder, filters.int.and-combination, filters.int.master-active-toggle, filters.int.active-counter-counts-filtering-only, filters.int.header-search-hides-cards, filters.int.esc-toggles-not-resets]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {cardCount, clickResetCriteriaIcon, driveOpenMenuLeaf, expectHeaderCounter, expectHeaderCounterNow,
  expectHeaderCounterQuiet, headerCounterTarget, trueCount} from '../../helpers/filter-panel';
import {addCardViaPicker} from './column-picker';
import {cardCaptions, categoryRowPoint, checkboxCensus, collapseHeaderSearch, driveHeaderSearch, filterState, FULL,
  headerSearchState, RACE_CATEGORY, removeAllViaPanelMenu, visibleCardCaptions, waitForPanelSettled} from './panel-core-ladder-shared';

declare const grok: any;

// Steps 1-9 of the ladder on the local lane. The layout and project round-trips (Steps 10, 11,
// 11-negative, 13) live in panel-core-ladder-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function toggleCardCheckbox(page: Page, column: string): Promise<{count: number; checked: boolean}> {
  return page.evaluate(async (column: string) => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === column);
    if (!card) throw new Error(`toggleCardCheckbox: the panel carries no ${column} card`);
    const cb = card.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null;
    if (!cb) throw new Error(`toggleCardCheckbox: the ${column} card carries no enable/disable checkbox`);
    const was = grok.shell.tv.dataFrame.filter.trueCount;
    const count = await (window as any).__filtered(() => cb.click(), 3000, was);
    return {count, checked: cb.checked};
  }, column);
}


async function tooltipStates(page: Page): Promise<Array<{
  display: string; visibility: string; text: string; cells: string[];
}>> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-tooltip')).map((t) => {
      const cs = window.getComputedStyle(t as HTMLElement);
      return {
        display: cs.display,
        visibility: cs.visibility,
        text: (t.textContent ?? '').trim(),
        cells: Array.from(t.querySelectorAll('td,th')).map((c) => c.textContent?.trim() ?? ''),
      };
    }));
}

async function raiseCounterTooltipCells(page: Page): Promise<string[]> {
  const target = await headerCounterTarget(page);
  if (!target.present)
    throw new Error('the header active-filter counter is not in the DOM — there is nothing to hover');
  if (!target.visible) {
    throw new Error('the header active-filter counter is not on screen to hover ' +
      `(display=${target.display}, visibility=${target.visibility}, box=${target.w}x${target.h}, ` +
      `text=${JSON.stringify(target.text)}) — the counter block is hidden while nothing is filtering`);
  }
  await page.mouse.move(5, 5);
  await v.pollValue(async () => (await tooltipStates(page))
    .filter((t) => t.display !== 'none' && t.visibility !== 'hidden').length,
  (shown) => shown === 0, 800, 50);
  const baseline = (await tooltipStates(page)).map((t) => t.text);
  const fresh = (i: number, text: string) => text.length > 0 && text !== (baseline[i] ?? '');
  await page.mouse.move(target.x, Math.max(2, target.y - 60), {steps: 6});
  await page.waitForTimeout(150);
  await page.mouse.move(target.x, target.y, {steps: 6});
  await page.waitForTimeout(150);
  await page.mouse.move(target.x, target.y + 1, {steps: 2});

  let states: Array<{display: string; visibility: string; text: string; cells: string[]}> = [];
  const shown = await v.pollValue(async () => {
    states = await tooltipStates(page);
    return states.findIndex((t, i) => t.display !== 'none' && t.visibility !== 'hidden' && fresh(i, t.text));
  }, (idx) => idx >= 0, 6000, 200);
  if (shown >= 0) return states[shown].cells;
  const after = await headerCounterTarget(page);
  const menus = await page.evaluate(() => document.querySelectorAll('.d4-menu-popup').length);
  const seen = states.length === 0 ? 'none' : states
    .map((t, i) => `#${i}: display=${t.display}, visibility=${t.visibility}, textLength=${t.text.length}, ` +
      `changed=${fresh(i, t.text)}`)
    .join('; ');
  throw new Error('the header counter tooltip did not come up within 6s of the hover. ' +
    `.d4-tooltip elements: ${seen}. Counter: present=${after.present}, visible=${after.visible}, ` +
    `display=${after.display}, visibility=${after.visibility}, text=${JSON.stringify(after.text)}, ` +
    `box=${after.w}x${after.h}. Hovered at (${Math.round(target.x)}, ${Math.round(target.y)}); ` +
    `open menu popups: ${menus}`);
}

test('Filter Panel — Panel Core Ladder', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  // the shared page can hand over a panel that reapplies an earlier spec's filter state
  await removeAllViaPanelMenu(page);
  await v.resetFilters(page);
  await expect.poll(() => trueCount(page), {timeout: 10_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(FULL);
  expect(await cardCaptions(page)).toEqual([]);
  expect(await cardCount(page)).toBe(0);

  await addCardViaPicker(page, 'SEX');
  expect(await cardCaptions(page)).toEqual(['SEX']);
  expect(await cardCount(page)).toBe(1);
  expect(await trueCount(page)).toBe(FULL);
  const counterAtSetup = await headerCounterTarget(page);
  expect(counterAtSetup.present).toBe(true);
  expect(counterAtSetup.visible,
    `the counter block is on screen while nothing filters (text ${JSON.stringify(counterAtSetup.text)})`).toBe(false);

  let trueCountRaceOnly = -1;
  let trueCountAgeOnly = -1;

  try {
    await softStep('Step 1 Drag RACE column header onto the panel → RACE card appears', async () => {
      const press = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const mainGrid = document.querySelector('.d4-table-view [name="viewer-Grid"]') ?? document.querySelector('[name="viewer-Grid"]');
        const overlay = grid.overlay ?? mainGrid!.querySelector('[name="overlay"]');
        const rect = overlay.getBoundingClientRect();
        const gcol = grid.columns.byName('RACE');
        return {x: rect.left + (gcol.left + gcol.right) / 2, y: rect.top + grid.colHeaderHeight / 2};
      });
      await page.mouse.move(press.x, press.y);
      await page.mouse.down();
      for (let i = 1; i <= 8; i++)
        await page.mouse.move(press.x + i * 18, press.y + i * 26, {steps: 1});
      const zone = await page.evaluate(() => {
        const z = Array.from(document.querySelectorAll('.d4-drop-zone'))
          .find((e) => e.textContent?.trim() === 'Add filter' && e.parentElement === document.body) as HTMLElement | undefined;
        if (!z) return null;
        const r = z.getBoundingClientRect();
        return {x: r.left + r.width / 2, y: r.top + r.height / 2};
      });
      if (zone) {
        await page.mouse.move(zone.x, zone.y, {steps: 3});
        await page.mouse.up();
      } else {
        await page.mouse.up();
      }
      const hasRace = await v.pollValue(async () => page.evaluate(() =>
        Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .some((c) => c.textContent?.trim() === 'RACE')), (seen) => seen, 900, 50);
      expect(hasRace).toBe(true);
      // the panel keeps its drop state a beat after the release, and a card in that state takes no click
      await v.pollValue(() => page.evaluate(() =>
        document.body.classList.contains('d4-drag') || document.querySelectorAll('.d4-drop-zone').length > 0),
      (dragging) => !dragging, 2000, 50);
      await page.mouse.move(0, 0);
    });

    await softStep('Step 3 Click a RACE category row on the card canvas → filter narrows to it, counter reads 1', async () => {
      const before = await trueCount(page);
      const pt = await categoryRowPoint(page, 'RACE', RACE_CATEGORY);
      expect(pt, 'the RACE card never painted a body tall enough to carry the category rows').not.toBeNull();
      // a real pointer gesture: the card resolves the click through the row under the pointer,
      // which a dispatched click on a freshly dropped card does not set
      await page.mouse.move(pt!.x, pt!.y, {steps: 4});
      await page.waitForTimeout(150);
      await page.mouse.click(pt!.x, pt!.y);
      const after = await v.pollValue(() => trueCount(page), (c) => c !== before, 3000, 50);
      expect(after, `the category click at (${Math.round(pt!.x)}, ${Math.round(pt!.y)}) left the rows at ${before}; ` +
        `RACE state ${JSON.stringify(await filterState(page, 'RACE', 'categorical'))}`).not.toBe(before);
      expect(after).toBeLessThan(FULL);
      expect(after).toBeGreaterThan(0);
      await expectHeaderCounter(page, '1', 'the category click left one filtering card, so the header counter must read 1');
      const raceAfterClick = await filterState(page, 'RACE', 'categorical');
      expect(raceAfterClick, 'the RACE card exposes no filter state after the category click')
        .not.toBeNull();
      expect(raceAfterClick!.selected,
        'the category click left no category selection on the RACE card').not.toBeNull();
      expect(raceAfterClick!.selected).toEqual([RACE_CATEGORY]);
      trueCountRaceOnly = after;
    });

    await softStep('Step 4 Add the AGE card via the panel header column selector, Step 5 window incl. ' +
      'out-of-range end → counter reads 2, slider still moves', async () => {
      await addCardViaPicker(page, 'AGE');

      const ageCard = page.locator('[name="viewer-Filters"] .d4-filter')
        .filter({has: page.locator('.d4-filter-column-name', {hasText: /^AGE$/})});
      const cardBox = await ageCard.first().boundingBox();
      expect(cardBox).not.toBeNull();
      await page.mouse.move(cardBox!.x + cardBox!.width / 2, cardBox!.y + 8, {steps: 6});
      const indicator = ageCard.locator('.d4-filter-indicator').first();
      await v.pollValue(() => indicator.isVisible(), (shown) => shown, 400, 50);
      if (await indicator.isVisible())
        await indicator.click();
      else
        await indicator.evaluate((el: HTMLElement) => el.click());
      await driveOpenMenuLeaf(page, null, 'Min / max');
      const ageMax = ageCard.locator('input.d4-filter-input-max').first();
      await ageMax.waitFor({state: 'visible', timeout: 10_000});

      const rowsAfter = (from: number) =>
        v.pollValue(() => trueCount(page), (c) => c !== from, 1200, 50);
      await ageMax.fill('60');
      await ageMax.press('Enter');
      const at60 = await rowsAfter(trueCountRaceOnly);
      expect(trueCountRaceOnly, 'Step 3 did not record the categorical-only row count').toBeGreaterThan(0);
      expect(at60).toBeLessThan(trueCountRaceOnly);
      expect(at60).toBeGreaterThan(0);
      await expectHeaderCounter(page, '2', 'RACE and AGE are both filtering, so the header counter must read 2');

      await ageMax.fill('999');
      await ageMax.press('Enter');
      const at999 = await rowsAfter(at60);
      await ageMax.fill('55');
      await ageMax.press('Enter');
      const atMoved = await rowsAfter(at999);
      await ageMax.blur();
      expect(atMoved).not.toBe(at999);
      expect(atMoved).toBeLessThan(trueCountRaceOnly);
      expect(atMoved).toBeGreaterThan(0);
      await expectHeaderCounter(page, '2', 'both cards are still filtering after the out-of-range entry, so the counter must still read 2');
    });

    await softStep('Step 6 Counter tooltip summary lists exactly RACE and AGE', async () => {
      const cells = await raiseCounterTooltipCells(page);
      expect(cells.length).toBeGreaterThan(0);
      const joined = cells.join(' | ');
      expect(joined).toContain('RACE');
      expect(joined).toContain('AGE');
      expect(joined).not.toContain('SEX');

      const criterionFor = (column: string): string => {
        const i = cells.findIndex((c) => c === column);
        if (i < 0 || i + 1 >= cells.length) {
          throw new Error(`Step 6: the summary carries no criterion cell next to ${column} — ` +
            `cells: ${JSON.stringify(cells)}`);
        }
        return cells[i + 1];
      };
      const unselectedCategories: string[] = await page.evaluate((clicked: string) =>
        (grok.shell.tv.dataFrame.col('RACE').categories as string[])
          .filter((c) => !!c && c !== clicked), RACE_CATEGORY);
      expect(unselectedCategories.length,
        'demog RACE exposes no second category to look for').toBeGreaterThan(0);
      const raceCriterion = criterionFor('RACE');
      expect(raceCriterion.length, 'the RACE row carries an empty criterion cell').toBeGreaterThan(0);
      expect(raceCriterion).toContain(RACE_CATEGORY);
      expect(unselectedCategories.filter((c) => raceCriterion.includes(c)),
        `the RACE criterion ${JSON.stringify(raceCriterion)} names categories the click left out`)
        .toEqual([]);
      const ageCriterion = criterionFor('AGE');
      expect(ageCriterion.length, 'the AGE row carries an empty criterion cell').toBeGreaterThan(0);
    });

    await softStep('Step 7 Disable RACE card checkbox → trueCount rises to AGE-only, counter reads 1', async () => {
      const {before, after, present, disabledClass, checked, active, filtering} = await page.evaluate(async () => {
        const w = window as any;
        const df = grok.shell.tv.dataFrame;
        const before = df.filter.trueCount;
        const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
        const race = cards.find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE')!;
        const cb = race.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement;
        const raceFilter = () => grok.shell.tv.getFiltersGroup().filters
          .find((f: any) => { try { return f.filterColumnName === 'RACE' && f.filterType === 'categorical'; } catch (_) { return false; } });
        cb.click();
        const active = await w.__poll(() => raceFilter()?.isActive, (a: any) => a === false, 3000, 50);
        const after = await w.__moved(() => df.filter.trueCount, before, 3000);
        const stillPresent = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .some((c) => c.textContent?.trim() === 'RACE');
        return {before, after, present: stillPresent, disabledClass: race.classList.contains('d4-filter-disabled'),
          checked: cb.checked, active, filtering: raceFilter()?.isFiltering};
      });
      trueCountAgeOnly = after;
      expect(active, `the RACE card's checkbox click did not switch its filter off (checkbox checked=${checked}, ` +
        `isActive=${active}, isFiltering=${filtering}, rows ${before} -> ${after})`).toBe(false);
      expect(after, `the RACE filter is off (isFiltering=${filtering}) but the rows did not come back`).toBeGreaterThan(before);
      expect(after).toBeLessThan(FULL);
      await expectHeaderCounter(page, '1', 'disabling the RACE card leaves AGE alone filtering, so the counter must drop to 1');
      expect(present).toBe(true);
      expect(disabledClass).toBe(true);

      const raceOff = await filterState(page, 'RACE', 'categorical');
      expect(raceOff, 'the disabled RACE card exposes no filter state to read its criterion from')
        .not.toBeNull();
      expect(raceOff!.selected,
        'the disabled RACE card kept no category selection — its criterion was erased, not kept')
        .not.toBeNull();
      expect(raceOff!.selected).toEqual([RACE_CATEGORY]);

      const back = await toggleCardCheckbox(page, 'RACE');
      expect(back.checked).toBe(true);
      expect(back.count, 'the re-enabled RACE card did not reproduce the two-filter row count')
        .toBe(before);
      await expectHeaderCounter(page, '2', 're-enabling the RACE card puts two cards back into filtering, so the counter must read 2');
      const offAgain = await toggleCardCheckbox(page, 'RACE');
      expect(offAgain.checked).toBe(false);
      expect(offAgain.count).toBe(trueCountAgeOnly);
      await expectHeaderCounter(page, '1', 'the RACE card is disabled again, so the counter must be back to 1');
    });

    await softStep('Step 8a Master toggle off → full count, d4-filters-disabled, cards present', async () => {
      const {after, disabled} = await page.evaluate(async () => {
        const panel = document.querySelector('[name="viewer-Filters"]')!;
        const was = grok.shell.tv.dataFrame.filter.trueCount;
        const after = await (window as any).__filtered(() =>
          (panel.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement).click(), 3000, was);
        return {after, disabled: panel.classList.contains('d4-filters-disabled')};
      });
      expect(after).toBe(FULL);
      expect(disabled).toBe(true);
      const captions = await cardCaptions(page);
      expect(captions).toContain('RACE');
      expect(captions).toContain('AGE');
      expect(captions).toContain('SEX');

      const raceState = await filterState(page, 'RACE', 'categorical');
      expect(raceState, 'the RACE card exposes no filter state while the master toggle is off')
        .not.toBeNull();
      expect(raceState!.selected,
        'the RACE card lost its category selection when the master toggle went off').not.toBeNull();
      expect(raceState!.selected).toEqual([RACE_CATEGORY]);
      const ageState = await filterState(page, 'AGE', 'histogram');
      expect(ageState, 'the AGE card exposes no filter state while the master toggle is off')
        .not.toBeNull();
      expect(ageState!.min, 'the AGE card lost its lower bound when the master toggle went off')
        .not.toBeNull();
      expect(ageState!.max, 'the AGE card lost its upper bound when the master toggle went off')
        .not.toBeNull();
      const ageColMax = await page.evaluate(() => grok.shell.tv.dataFrame.col('AGE').max);
      expect(ageState!.max!).toBeGreaterThan(ageState!.min!);
      expect(ageState!.max!).toBeLessThan(ageColMax);
    });

    await softStep('Step 8b Master toggle on → restores AGE-only value, RACE checkbox stays unchecked', async () => {
      const {after, disabled, raceChecked} = await page.evaluate(async () => {
        const panel = document.querySelector('[name="viewer-Filters"]')!;
        const was = grok.shell.tv.dataFrame.filter.trueCount;
        const after = await (window as any).__filtered(() =>
          (panel.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement).click(), 3000, was);
        const race = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
          .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE')!;
        return {
          after,
          disabled: panel.classList.contains('d4-filters-disabled'),
          raceChecked: (race.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement).checked,
        };
      });
      expect(after).toBe(trueCountAgeOnly);
      expect(disabled).toBe(false);
      expect(raceChecked).toBe(false);
    });

    await softStep('Step 8c Esc toggles the group off and back on, criteria intact throughout', async () => {
      const before = await trueCount(page);
      expect(before, 'Esc has nothing to toggle unless the panel is filtering on the way in')
        .toBeLessThan(FULL);
      expect(before).toBe(trueCountAgeOnly);

      const panel = page.locator('[name="viewer-Filters"]').first();
      await panel.click({position: {x: 5, y: 5}});
      await page.keyboard.press('Escape');
      await expect.poll(async () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
        {timeout: 10_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(FULL);

      const off = await page.evaluate(() => ({
        rows: grok.shell.tv.dataFrame.filter.trueCount,
        disabled: document.querySelector('[name="viewer-Filters"]')!.classList.contains('d4-filters-disabled'),
      }));
      expect(off.rows).toBe(FULL);
      expect(off.disabled, 'Esc did not put the group into its disabled state').toBe(true);
      const raceAfterEsc = await filterState(page, 'RACE', 'categorical');
      expect(raceAfterEsc, 'the RACE card lost its state to Esc — that is reset behaviour, not toggle')
        .not.toBeNull();
      expect(raceAfterEsc!.selected).toEqual([RACE_CATEGORY]);
      const ageAfterEsc = await filterState(page, 'AGE', 'histogram');
      expect(ageAfterEsc, 'the AGE card lost its state to Esc').not.toBeNull();
      expect(ageAfterEsc!.max, 'Esc cleared the AGE window — a reset would, a toggle must not')
        .not.toBeNull();

      await panel.click({position: {x: 5, y: 5}});
      await page.keyboard.press('Escape');
      await expect.poll(async () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
        {timeout: 10_000, intervals: [30, 60, 120, 250, 500, 1000]}).toBe(before);
      expect(await page.evaluate(() =>
        document.querySelector('[name="viewer-Filters"]')!.classList.contains('d4-filters-disabled')))
        .toBe(false);
      await expectHeaderCounter(page, '1', 'the counter did not come back with the filtering after the second Esc');
    });

    await softStep('Step 9 Reset icon → full count, counter 0, cards present, checkboxes re-checked', async () => {
      const counterBeforeSearch = await headerCounterTarget(page);
      expect(counterBeforeSearch.present,
        'the header active-filter counter is not in the DOM, so there is nothing to compare across ' +
        'the search — the check below would be vacuous').toBe(true);
      expect(counterBeforeSearch.visible,
        'the header active-filter counter is off screen while a card is filtering, so the ' +
        'comparison across the search would be vacuous').toBe(true);
      expect(counterBeforeSearch.text,
        'the counter is on screen but blank, so the comparison across the search would be vacuous')
        .not.toBe('');
      const searched = await driveHeaderSearch(page, 'RACE');
      expect(searched.typed).toBe('RACE');
      expect(searched.visibleAfter.length).toBeLessThan(searched.visibleBefore.length);
      expect(searched.visibleAfter).toContain('RACE');
      expect(searched.visibleAfter.every((c) => c.toUpperCase().includes('RACE'))).toBe(true);
      expect(searched.visibleAfter,
        'the AGE card — the one that is actually filtering — was not hidden, so the counter check ' +
        'below is not about a hidden filtering card').not.toContain('AGE');
      expect(searched.countAfter).toBe(searched.countBefore);
      await expectHeaderCounterNow(page, counterBeforeSearch.text,
        'the active-filter counter moved when a filtering card was merely hidden by the header search');
      await collapseHeaderSearch(page);
      expect((await visibleCardCaptions(page)).length).toBe(searched.visibleBefore.length);
      expect(await trueCount(page)).toBe(searched.countBefore);

      expect(searched.countBefore,
        'nothing is filtering on the way into the reset, so the settle barrier below could not tell ' +
        'a completed reset from a stale reading').toBeLessThan(FULL);
      await clickResetCriteriaIcon(page, {via: 'dom'});
      const rows = await waitForPanelSettled(page, {changedFrom: searched.countBefore, timeoutMs: 10_000});
      expect(rows).toBe(FULL);
      expect(await page.evaluate(() => !!document.querySelector('.d4-dialog'))).toBe(false);
      await expectHeaderCounterQuiet(page,
        'after the header reset the counter must be hidden or read 0');
      const captions = await cardCaptions(page);
      expect(captions).toContain('RACE');
      expect(captions).toContain('AGE');
      expect(captions).toContain('SEX');
      const census = await checkboxCensus(page);
      expect(census.boxes).toBeGreaterThan(0);
      expect(census.boxes).toBe(census.cards);
      expect(census.checked).toBe(census.boxes);
      expect((await visibleCardCaptions(page)).length).toBe(census.cards);
      const search = await headerSearchState(page);
      if (search === null)
        throw new Error('Step 9: the panel header search input is not in the DOM — its state cannot be read');
      expect(search.visible ? search.value : '',
        `header search after reset: visible=${search.visible}, value=${JSON.stringify(search.value)}`).toBe('');
    });
  } finally {
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
