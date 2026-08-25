/* ---
realizes: [filters.cp.panel-core-ladder, filters.int.and-combination, filters.int.master-active-toggle, filters.int.active-counter-counts-filtering-only, filters.int.header-search-hides-cards, filters.int.esc-toggles-not-resets]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import {addCardViaColumnSelector, cardCount, clickResetCriteriaIcon, driveOpenMenuLeaf,
  expectHeaderCounter, expectHeaderCounterNow, expectHeaderCounterQuiet, headerCounterTarget,
  trueCount} from '../../helpers/filter-panel';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const FULL = 5850;
const RACE_CATEGORY = 'Black';
const projectName = `filters-ladder-${Date.now()}`;
const runTag = `ladder-${Date.now()}-${Math.floor(Math.random() * 1e9)}`;

async function cardCaptions(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
      .map((c) => c.textContent?.trim() ?? ''));
}

async function visibleCardCaptions(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .filter((card) => (card as HTMLElement).offsetParent !== null)
      .map((card) => card.querySelector('.d4-filter-column-name')?.textContent?.trim() ?? ''));
}

async function checkboxCensus(page: import('@playwright/test').Page):
    Promise<{cards: number; boxes: number; checked: number}> {
  return page.evaluate(() => {
    const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
    const boxes = cards
      .map((c) => c.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null)
      .filter((cb) => cb !== null) as HTMLInputElement[];
    return {cards: cards.length, boxes: boxes.length, checked: boxes.filter((cb) => cb.checked).length};
  });
}

async function filterState(page: import('@playwright/test').Page, column: string, type: string):
    Promise<{selected: string[] | null; min: number | null; max: number | null} | null> {
  return page.evaluate(({column, type}) => {
    const states = grok.shell.tv.getFiltersGroup().getStates(column, type) as any[];
    if (!states || states.length === 0) return null;
    const s = states[0];
    return {
      selected: Array.isArray(s.selected) ? s.selected.map((c: any) => String(c)) : null,
      min: typeof s.min === 'number' ? s.min : null,
      max: typeof s.max === 'number' ? s.max : null,
    };
  }, {column, type});
}

async function toggleCardCheckbox(page: import('@playwright/test').Page, column: string):
    Promise<{count: number; checked: boolean}> {
  return page.evaluate(async (column: string) => {
    const card = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === column);
    if (!card) throw new Error(`toggleCardCheckbox: the panel carries no ${column} card`);
    const cb = card.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null;
    if (!cb) throw new Error(`toggleCardCheckbox: the ${column} card carries no enable/disable checkbox`);
    cb.click();
    await new Promise((r) => setTimeout(r, 800));
    return {count: grok.shell.tv.dataFrame.filter.trueCount, checked: cb.checked};
  }, column);
}

async function waitForPanelSettled(page: import('@playwright/test').Page,
  opts: {changedFrom?: number; timeoutMs?: number} = {}): Promise<number> {
  const timeoutMs = opts.timeoutMs ?? 60_000;
  const deadline = Date.now() + timeoutMs;
  let stable = 0;
  let last = Number.NaN;
  for (;;) {
    const now = await page.evaluate(() => ({
      ready: !!document.querySelector('[name="viewer-Filters"]') &&
        document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length > 0,
      count: grok.shell.tv?.dataFrame?.filter?.trueCount ?? -1,
    }));
    const moved = opts.changedFrom === undefined || now.count !== opts.changedFrom;
    stable = now.ready && now.count === last ? stable + 1 : 0;
    last = now.count;
    if (now.ready && moved && stable >= 2) return now.count;
    if (Date.now() > deadline) {
      throw new Error(`waitForPanelSettled: the panel never settled within ${timeoutMs}ms ` +
        `(last row count ${last}, waiting for it to leave ${opts.changedFrom})`);
    }
    await page.waitForTimeout(500);
  }
}

async function driveHeaderSearch(page: import('@playwright/test').Page, text: string): Promise<{
  visibleBefore: string[]; visibleAfter: string[]; countBefore: number; countAfter: number; typed: string;
}> {
  const visibleBefore = await visibleCardCaptions(page);
  const countBefore = await trueCount(page);
  const input = page.locator('[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]');
  if (!(await input.isVisible())) {
    const opened = await page.evaluate(() => {
      const icon = document.querySelector(
        '[name="viewer-Filters"] .d4-filter-group-header [name="icon-search"]') as HTMLElement | null;
      if (!icon) return false;
      icon.click();
      return true;
    });
    if (!opened)
      throw new Error('driveHeaderSearch: the panel header carries no search icon — the search cannot be opened');
  }
  await input.waitFor({state: 'visible', timeout: 15_000});
  await input.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.type(text, {delay: 40});
  const typed = await input.inputValue();
  await page.waitForTimeout(800);
  return {
    visibleBefore,
    visibleAfter: await visibleCardCaptions(page),
    countBefore,
    countAfter: await trueCount(page),
    typed,
  };
}

async function collapseHeaderSearch(page: import('@playwright/test').Page): Promise<void> {
  const input = page.locator('[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]');
  if (!(await input.isVisible()))
    throw new Error('collapseHeaderSearch: the header search box is not open, so it cannot be closed');
  await page.evaluate(() => {
    const icon = document.querySelector(
      '[name="viewer-Filters"] .d4-filter-group-header [name="icon-search"]') as HTMLElement | null;
    if (!icon) throw new Error('collapseHeaderSearch: the panel header carries no search icon');
    icon.click();
  });
  await input.waitFor({state: 'hidden', timeout: 10_000});
  await page.waitForTimeout(400);
}

async function headerSearchState(page: import('@playwright/test').Page):
    Promise<{value: string; visible: boolean} | null> {
  return page.evaluate(() => {
    const el = document.querySelector(
      '[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]') as HTMLInputElement | null;
    return el === null ? null : {value: el.value, visible: el.offsetParent !== null};
  });
}

async function removeAllViaPanelMenu(page: import('@playwright/test').Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(async () => cardCount(page),
    {timeout: 20_000, intervals: [300, 600, 1200],
      message: 'the panel menu\'s "Remove All" leaf was driven but the panel still carries cards'})
    .toBe(0);
}

async function myApplicableLayouts(page: import('@playwright/test').Page):
    Promise<Array<{id: string; name: string}>> {
  return page.evaluate(async () => {
    const me = String(grok.shell.user.id);
    const ls = (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)) ?? [];
    return ls
      .filter((l: any) => !l.author || !l.author.id || String(l.author.id) === me)
      .map((l: any) => ({id: String(l.id), name: String(l.friendlyName ?? l.name ?? '')}));
  });
}

function carriesMarker(name: string, marker: string): boolean {
  return name === marker || name.startsWith(`${marker} (`);
}

async function stampTableName(page: import('@playwright/test').Page, marker: string): Promise<void> {
  const applied = await page.evaluate((name: string) => {
    grok.shell.tv.dataFrame.name = name;
    return String(grok.shell.tv.dataFrame.name);
  }, marker);
  expect(applied, 'the open table did not take this run\'s marker name, so a layout saved from it ' +
    'could not be told apart from one another process saved as this same user').toBe(marker);
}

async function saveLayoutToGallery(page: import('@playwright/test').Page, marker: string): Promise<string> {
  const stamped = await page.evaluate(() => String(grok.shell.tv?.dataFrame?.name ?? ''));
  expect(stamped, 'the table is not carrying this run\'s marker at save time, so the layout the ' +
    'save produces would be unattributable').toBe(marker);

  const beforeIds = (await myApplicableLayouts(page)).map((l) => l.id);
  expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery']),
    'the View | Layout | Save to Gallery leaf could not be driven, so no layout was saved').toBe(true);

  let mine: Array<{id: string; name: string}> = [];
  await expect.poll(async () => {
    mine = (await myApplicableLayouts(page))
      .filter((l) => !beforeIds.includes(l.id) && carriesMarker(l.name, marker));
    return mine.length;
  }, {
    timeout: 25_000,
    intervals: [500, 1000, 2000, 3000],
    message: `the gallery save produced no new applicable layout named after this run's marker "${marker}"`,
  }).toBeGreaterThanOrEqual(1);
  if (mine.length !== 1) {
    throw new Error(
      `The gallery save produced ${mine.length} new applicable layouts carrying this run's marker ` +
      `"${marker}" (${mine.map((l) => `${l.id}:${l.name}`).join(', ')}), expected exactly 1 — ` +
      'refusing to guess which one is ours; deleting none.');
  }
  return mine[0].id;
}

async function categoryRowPoint(page: import('@playwright/test').Page, column: string, category: string) {
  return page.evaluate(({column, category}) => {
    const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
    const card = cards.find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === column);
    if (!card) return null;
    const canvas = card.querySelector('[name="canvas"]') as HTMLCanvasElement | null;
    if (!canvas) return null;
    const cats: string[] = grok.shell.tv.dataFrame.col(column).categories;
    const idx = cats.indexOf(category);
    if (idx < 0) return null;
    const rect = canvas.getBoundingClientRect();
    const rowHeight = rect.height / cats.length;
    return {x: rect.left + 12, y: rect.top + idx * rowHeight + rowHeight / 2};
  }, {column, category});
}

async function tooltipStates(page: import('@playwright/test').Page): Promise<Array<{
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

async function raiseCounterTooltipCells(page: import('@playwright/test').Page): Promise<string[]> {
  const target = await headerCounterTarget(page);
  if (!target.present)
    throw new Error('the header active-filter counter is not in the DOM — there is nothing to hover');
  if (!target.visible) {
    throw new Error('the header active-filter counter is not on screen to hover ' +
      `(display=${target.display}, visibility=${target.visibility}, box=${target.w}x${target.h}, ` +
      `text=${JSON.stringify(target.text)}) — the counter block is hidden while nothing is filtering`);
  }
  await page.mouse.move(5, 5);
  await page.waitForTimeout(800);
  const baseline = (await tooltipStates(page)).map((t) => t.text);
  const fresh = (i: number, text: string) => text.length > 0 && text !== (baseline[i] ?? '');
  await page.mouse.move(target.x, Math.max(2, target.y - 60), {steps: 6});
  await page.waitForTimeout(150);
  await page.mouse.move(target.x, target.y, {steps: 6});
  await page.waitForTimeout(150);
  await page.mouse.move(target.x, target.y + 1, {steps: 2});

  const deadline = Date.now() + 6_000;
  let states: Array<{display: string; visibility: string; text: string; cells: string[]}> = [];
  for (;;) {
    states = await tooltipStates(page);
    const idx = states.findIndex((t, i) =>
      t.display !== 'none' && t.visibility !== 'hidden' && fresh(i, t.text));
    if (idx >= 0) return states[idx].cells;
    if (Date.now() > deadline) break;
    await page.waitForTimeout(200);
  }
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

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  expect(await trueCount(page)).toBe(FULL);
  await removeAllViaPanelMenu(page);
  expect(await cardCaptions(page)).toEqual([]);
  expect(await cardCount(page)).toBe(0);
  expect(await trueCount(page)).toBe(FULL);

  await addCardViaColumnSelector(page, 'SEX');
  expect(await cardCaptions(page)).toEqual(['SEX']);
  expect(await cardCount(page)).toBe(1);
  expect(await trueCount(page)).toBe(FULL);
  const counterAtSetup = await headerCounterTarget(page);
  expect(counterAtSetup.present).toBe(true);
  expect(counterAtSetup.visible,
    `the counter block is on screen while nothing filters (text ${JSON.stringify(counterAtSetup.text)})`).toBe(false);

  let trueCountRaceOnly = -1;
  let trueCountAgeOnly = -1;
  let trueCountSaved = -1;
  let projectId = '';

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
      await page.waitForTimeout(900);
      const hasRace = await page.evaluate(() =>
        Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .some((c) => c.textContent?.trim() === 'RACE'));
      expect(hasRace).toBe(true);
    });

    await softStep('Step 3 Click a RACE category row on the card canvas → filter narrows to it, counter reads 1', async () => {
      const before = await trueCount(page);
      const pt = await categoryRowPoint(page, 'RACE', RACE_CATEGORY);
      expect(pt).not.toBeNull();
      await page.mouse.click(pt!.x, pt!.y);
      await page.waitForTimeout(900);
      const after = await trueCount(page);
      expect(after).not.toBe(before);
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
      await addCardViaColumnSelector(page, 'AGE');

      const ageCard = page.locator('[name="viewer-Filters"] .d4-filter')
        .filter({has: page.locator('.d4-filter-column-name', {hasText: /^AGE$/})});
      const cardBox = await ageCard.first().boundingBox();
      expect(cardBox).not.toBeNull();
      await page.mouse.move(cardBox!.x + cardBox!.width / 2, cardBox!.y + 8, {steps: 6});
      await page.waitForTimeout(400);
      const indicator = ageCard.locator('.d4-filter-indicator').first();
      if (await indicator.isVisible())
        await indicator.click();
      else
        await indicator.evaluate((el: HTMLElement) => el.click());
      await driveOpenMenuLeaf(page, null, 'Min / max');
      const ageMax = ageCard.locator('input.d4-filter-input-max').first();
      await ageMax.waitFor({state: 'visible', timeout: 10_000});

      const rowsNow = async () => trueCount(page);
      await ageMax.fill('60');
      await ageMax.press('Enter');
      await page.waitForTimeout(1200);
      const at60 = await rowsNow();
      expect(trueCountRaceOnly, 'Step 3 did not record the categorical-only row count').toBeGreaterThan(0);
      expect(at60).toBeLessThan(trueCountRaceOnly);
      expect(at60).toBeGreaterThan(0);
      await expectHeaderCounter(page, '2', 'RACE and AGE are both filtering, so the header counter must read 2');

      await ageMax.fill('999');
      await ageMax.press('Enter');
      await page.waitForTimeout(1200);
      const at999 = await rowsNow();
      await ageMax.fill('55');
      await ageMax.press('Enter');
      await page.waitForTimeout(1200);
      const atMoved = await rowsNow();
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
      const {before, after, present, disabledClass} = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const before = df.filter.trueCount;
        const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
        const race = cards.find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE')!;
        (race.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement).click();
        await new Promise((r) => setTimeout(r, 800));
        const stillPresent = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .some((c) => c.textContent?.trim() === 'RACE');
        return {before, after: df.filter.trueCount, present: stillPresent, disabledClass: race.classList.contains('d4-filter-disabled')};
      });
      trueCountAgeOnly = after;
      expect(after).toBeGreaterThan(before);
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
        const df = grok.shell.tv.dataFrame;
        const panel = document.querySelector('[name="viewer-Filters"]')!;
        (panel.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement).click();
        await new Promise((r) => setTimeout(r, 800));
        return {
          after: df.filter.trueCount,
          disabled: panel.classList.contains('d4-filters-disabled'),
        };
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
        const df = grok.shell.tv.dataFrame;
        const panel = document.querySelector('[name="viewer-Filters"]')!;
        (panel.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement).click();
        await new Promise((r) => setTimeout(r, 800));
        const race = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
          .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE')!;
        return {
          after: df.filter.trueCount,
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
        {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(FULL);

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
        {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(before);
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
      const rows = await waitForPanelSettled(page, {changedFrom: searched.countBefore});
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

    await softStep('Step 10 Re-establish the two-filter state and save the layout, Step 11 layout ' +
      'round-trip → filters restored to the saved value', async () => {
      const mainMarker = `${runTag}-main`;
      let savedLayout: any = null;
      try {
        await stampTableName(page, mainMarker);
        trueCountSaved = await page.evaluate(async (category: string) => {
          const fg = grok.shell.tv.getFiltersGroup();
          const df = grok.shell.tv.dataFrame;
          fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: [category]});
          await new Promise((r) => setTimeout(r, 600));
          fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 18, max: 60});
          await new Promise((r) => setTimeout(r, 700));
          return df.filter.trueCount;
        }, RACE_CATEGORY);
        await expectHeaderCounter(page, '2', 'the two-filter state is re-established, so the counter must read 2 before the layout is saved');

        savedLayout = await saveLayoutToGallery(page, mainMarker);

        await page.evaluate(async () => {
          grok.shell.tv.dataFrame.resetFilter();
          try { grok.shell.tv.addViewer('Bar chart'); } catch (_) {}
          await new Promise((r) => setTimeout(r, 800));
        });
        const perturbed = await trueCount(page);
        expect(perturbed).not.toBe(trueCountSaved);

        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          grok.shell.tv.loadLayout(saved);
        }, savedLayout);
        const settled = await waitForPanelSettled(page, {changedFrom: perturbed});
        const caps = await cardCaptions(page);
        expect(await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'))).toBe(true);
        expect(caps).toContain('RACE');
        expect(caps).toContain('AGE');
        const census = await checkboxCensus(page);
        expect(census.boxes).toBeGreaterThan(0);
        expect(census.boxes).toBe(census.cards);
        expect(census.checked).toBe(census.boxes);
        const raceRestored = await filterState(page, 'RACE', 'categorical');
        expect(raceRestored, 'the re-applied layout left the RACE card with no filter state')
          .not.toBeNull();
        expect(raceRestored!.selected,
          'the re-applied RACE card carries no category selection — the criterion did not survive')
          .not.toBeNull();
        expect(raceRestored!.selected).toEqual([RACE_CATEGORY]);
        const ageRestored = await filterState(page, 'AGE', 'histogram');
        expect(ageRestored, 'the re-applied layout left the AGE card with no filter state')
          .not.toBeNull();
        expect(ageRestored!.max, 'the re-applied AGE card carries no upper bound — its window was lost')
          .not.toBeNull();
        const ageColMaxAfterLayout = await page.evaluate(() => grok.shell.tv.dataFrame.col('AGE').max);
        expect(ageRestored!.max!).toBeLessThan(ageColMaxAfterLayout);
        expect(settled).toBe(trueCountSaved);
      } finally {
        if (savedLayout) {
          await page.evaluate(async (layoutId: string) => {
            try { const s = await grok.dapi.layouts.find(layoutId); await grok.dapi.layouts.delete(s); } catch (_) {}
          }, savedLayout);
        }
      }
    });

    await softStep('Step 11-negative GROK-16677 reset-then-save layout → full count, panel search empty', async () => {
      const resetMarker = `${runTag}-reset`;
      let resetLayout: any = null;
      try {
        await stampTableName(page, resetMarker);
        const cardsBefore = await cardCount(page);
        const countBeforeReset = await trueCount(page);
        expect(countBeforeReset,
          'nothing is filtering on the way into this reset, so the settle barrier below could not ' +
          'tell a completed reset from a stale reading').toBeLessThan(FULL);
        await clickResetCriteriaIcon(page, {via: 'dom'});
        const resetCount = await waitForPanelSettled(page, {changedFrom: countBeforeReset});
        expect(cardsBefore).toBeGreaterThan(0);
        expect(await cardCount(page)).toBe(cardsBefore);
        expect(resetCount).toBe(FULL);

        const searched = await driveHeaderSearch(page, 'RACE');
        expect(searched.typed).toBe('RACE');
        expect(searched.visibleAfter.length).toBeLessThan(searched.visibleBefore.length);
        expect(searched.countAfter).toBe(searched.countBefore);

        resetLayout = await saveLayoutToGallery(page, resetMarker);

        const perturbed = await page.evaluate(async () => {
          const df = grok.shell.tv.dataFrame;
          grok.shell.tv.getFiltersGroup()
            .updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
          await new Promise((r) => setTimeout(r, 900));
          return df.filter.trueCount;
        });
        expect(perturbed).not.toBe(FULL);
        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          grok.shell.tv.loadLayout(saved);
        }, resetLayout);
        const settled = await waitForPanelSettled(page, {changedFrom: perturbed});
        expect(settled).toBe(FULL);
        const search = await headerSearchState(page);
        if (search === null) {
          throw new Error(
            'Step 11-negative: the re-applied panel has no header search input — its value cannot be read');
        }
        expect(search.value).toBe('');
        const census = await checkboxCensus(page);
        expect(census.cards).toBeGreaterThan(0);
        expect((await visibleCardCaptions(page)).length).toBe(census.cards);
      } finally {
        if (resetLayout) {
          await page.evaluate(async (layoutId: string) => {
            try { const s = await grok.dapi.layouts.find(layoutId); await grok.dapi.layouts.delete(s); } catch (_) {}
          }, resetLayout);
        }
      }
    });

    await softStep('Step 13 Project round-trip → panel reopens, filters restored (GROK-19152 barrier)', async () => {
      trueCountSaved = await page.evaluate(async (category: string) => {
        const fg = grok.shell.tv.getFiltersGroup();
        const df = grok.shell.tv.dataFrame;
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: [category]});
        await new Promise((r) => setTimeout(r, 600));
        fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 18, max: 60});
        await new Promise((r) => setTimeout(r, 700));
        return df.filter.trueCount;
      }, RACE_CATEGORY);

      const saved = await saveProjectViaApi(page, projectName);
      projectId = saved.projectId;

      await page.evaluate(async () => { grok.shell.closeAll(); await new Promise((r) => setTimeout(r, 800)); });

      await page.evaluate(async (id: string) => {
        const project = await grok.dapi.projects.find(id);
        await project.open();
      }, projectId);
      await page.waitForFunction(() => {
        const el = document.querySelector('[name="viewer-Filters"]') as HTMLElement | null;
        return !!el && el.offsetParent !== null;
      }, {timeout: 60_000});

      const reopened = await page.evaluate(async () => {
        await new Promise((r) => setTimeout(r, 1500));
        const caps = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .map((c) => c.textContent?.trim());
        return {
          panel: !!document.querySelector('[name="viewer-Filters"]'),
          hasRace: caps.includes('RACE'),
          hasAge: caps.includes('AGE'),
          count: grok.shell.tv.dataFrame.filter.trueCount,
        };
      });
      expect(reopened.panel).toBe(true);
      expect(reopened.hasRace).toBe(true);
      expect(reopened.hasAge).toBe(true);
      const raceReopened = await filterState(page, 'RACE', 'categorical');
      expect(raceReopened, 'the reopened project left the RACE card with no filter state')
        .not.toBeNull();
      expect(raceReopened!.selected,
        'the reopened RACE card carries no category selection — the criterion was not persisted')
        .not.toBeNull();
      expect(raceReopened!.selected).toEqual([RACE_CATEGORY]);
      const ageReopened = await filterState(page, 'AGE', 'histogram');
      expect(ageReopened, 'the reopened project left the AGE card with no filter state').not.toBeNull();
      expect(ageReopened!.max, 'the reopened AGE card carries no upper bound — its window was not persisted')
        .not.toBeNull();
      const ageColMaxReopened = await page.evaluate(() => grok.shell.tv.dataFrame.col('AGE').max);
      expect(ageReopened!.max!).toBeLessThan(ageColMaxReopened);
      expect(reopened.count).toBe(trueCountSaved);
    });
  } finally {
    await deleteProjectWithCleanup(page, {projectId});
  }

  v.finishSpec();
});
