/* ---
realizes: [filters.cp.cloned-view-sync, filters.int.same-column-sync]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const FULL = 5850;

async function installViewResolver(page: import('@playwright/test').Page): Promise<void> {
  await page.evaluate(() => {
    (window as any).__tv = (vn: string) => {
      const tvs: any[] = [];
      for (const x of grok.shell.tableViews) tvs.push(x);
      const hits = tvs.filter((t: any) => t.name === vn);
      if (hits.length !== 1)
        throw new Error(`TableView "${vn}" matched ${hits.length} views (open: ${tvs.map((t: any) => t.name).join(' | ')})`);
      return hits[0];
    };
    (window as any).__cards = (vn: string): HTMLElement[] => {
      const root = (window as any).__tv(vn).root as HTMLElement;
      return Array.from(root.querySelectorAll('[name="viewer-Filters"] .d4-filter')) as HTMLElement[];
    };
    (window as any).__card = (vn: string, caption: string): HTMLElement => {
      const hits = (window as any).__cards(vn).filter((x: HTMLElement) =>
        x.querySelector('.d4-filter-column-name')?.textContent?.trim() === caption);
      if (hits.length !== 1)
        throw new Error(`the Filter Panel of view "${vn}" paints ${hits.length} cards captioned "${caption}", expected exactly one`);
      return hits[0];
    };
  });
}

async function viewNames(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() => {
    const names: string[] = [];
    for (const x of grok.shell.tableViews) names.push(x.name);
    return names;
  });
}

async function readFilters(page: import('@playwright/test').Page, viewName: string): Promise<
  {col: string | null; type: string; isActive: boolean; isFiltering: boolean}[]> {
  return page.evaluate((vn: string) => {
    return (window as any).__tv(vn).getFiltersGroup().filters.map((f: any) => {
      let col: string | null = null;
      try { col = f.filterColumnName ?? null; } catch (_) { /* proxy */ }
      if (!col) try { col = f.columnName ?? null; } catch (_) { /* proxy */ }
      return {col, type: f.filterType, isActive: f.isActive, isFiltering: f.isFiltering};
    });
  }, viewName);
}

async function panelCaptionsOf(page: import('@playwright/test').Page, viewName: string): Promise<string[]> {
  await activateView(page, viewName);
  return page.evaluate((vn: string) => {
    const root = (window as any).__tv(vn).root as HTMLElement;
    if (root.querySelector('[name="viewer-Filters"]') == null)
      throw new Error(`view "${vn}" has no Filter Panel in its own root`);
    return Array.from(root.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
      .map((e) => (e.textContent ?? '').trim());
  }, viewName);
}

async function cardDisabledIn(
  page: import('@playwright/test').Page, viewName: string, column: string,
): Promise<boolean> {
  await activateView(page, viewName);
  return page.evaluate(({vn, col}) =>
    (window as any).__card(vn, col).classList.contains('d4-filter-disabled'),
  {vn: viewName, col: column});
}

async function clickCardCheckboxIn(
  page: import('@playwright/test').Page, viewName: string, column: string,
): Promise<void> {
  await activateView(page, viewName);
  await page.evaluate(({vn, col}) => {
    const box = (window as any).__card(vn, col)
      .querySelector('input[type="checkbox"].ui-input-editor') as HTMLElement | null;
    if (!box)
      throw new Error(`the "${col}" card in view "${vn}" carries no enable/disable checkbox — the toggle gesture never happened`);
    box.click();
  }, {vn: viewName, col: column});
}

async function trueCountOf(page: import('@playwright/test').Page, viewName: string): Promise<number> {
  return page.evaluate((vn: string) => (window as any).__tv(vn).dataFrame.filter.trueCount, viewName);
}

async function selectedOf(
  page: import('@playwright/test').Page, viewName: string, column: string,
): Promise<string[]> {
  return page.evaluate(({vn, col}) => {
    const s = (window as any).__tv(vn).getFiltersGroup().getStates(col, 'categorical') as any[];
    if (!s || !s[0] || !Array.isArray(s[0].selected))
      throw new Error(`view "${vn}" has no categorical filter state for column "${col}"`);
    return s[0].selected.map((c: any) => String(c));
  }, {vn: viewName, col: column});
}

const DOT = 'icon-dot-circle';
const CIRCLE = 'icon-circle';
const MV_LEAVES: [string, string][] = [
  ['keep', 'div-Missing-values---Keep-missing-values'],
  ['filterOut', 'div-Missing-values---Filter-out-missing-values'],
  ['showOnly', 'div-Missing-values---Show-only-missing-value'],
];

function ageCard(page: import('@playwright/test').Page) {
  return page.locator('[name="viewer-Filters"]:visible .d4-filter')
    .filter({has: page.locator('.d4-filter-column-name', {hasText: /^AGE$/})});
}

async function boxOf(
  page: import('@playwright/test').Page, selector: string,
): Promise<{x: number; y: number; width: number; height: number}> {
  const box = await page.locator(selector).first().boundingBox();
  if (!box)
    throw new Error(`"${selector}" has no layout box — it is absent from the page or collapsed to zero size`);
  return box;
}

async function panelCounts(
  page: import('@playwright/test').Page, viewName: string,
): Promise<{cards: number; filters: number; ageCards: number}> {
  return page.evaluate((vn: string) => {
    const view = (window as any).__tv(vn);
    const cards = (window as any).__cards(vn) as HTMLElement[];
    return {
      cards: cards.length,
      filters: view.getFiltersGroup().filters.length,
      ageCards: cards.filter((c) =>
        c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'AGE').length,
    };
  }, viewName);
}

const CAT_ROW_TOP = 10;
const CAT_ROW_PITCH = 27;
const CAT_ROW_CENTRE = 13;
const CAT_NAME_X = 60;

async function activateView(page: import('@playwright/test').Page, viewName: string): Promise<void> {
  await page.evaluate((vn: string) => { grok.shell.v = (window as any).__tv(vn); }, viewName);
  await expect.poll(async () => page.evaluate((vn: string) => {
    const el = ((window as any).__tv(vn).root as HTMLElement)
      .querySelector('[name="viewer-Filters"]') as HTMLElement | null;
    return el != null && el.offsetParent !== null;
  }, viewName), {
    message: `the Filter Panel of view "${viewName}" never became the laid-out one after switching to it — `
      + 'every :visible panel read below would have addressed the other view',
    timeout: 15_000, intervals: [300, 600, 1200],
  }).toBe(true);
}

async function clickCategoryNameRowIn(
  page: import('@playwright/test').Page, viewName: string, column: string, rowIndex: number,
): Promise<void> {
  await activateView(page, viewName);
  await page.evaluate(({vn, col, x, y}) => {
    const overlay = (window as any).__card(vn, col)
      .querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null;
    if (!overlay)
      throw new Error(`the "${col}" card in view "${vn}" exposes no categorical [name="viewer-Grid"] [name="overlay"] body — the category-row click never happened`);
    const rect = overlay.getBoundingClientRect();
    if (rect.height <= y)
      throw new Error(`the "${col}" card body in view "${vn}" is ${rect.height}px tall, so the row at y=${y} is not painted and the click would land on empty canvas`);
    const o = {bubbles: true, cancelable: true, view: window, button: 0,
      clientX: rect.left + x, clientY: rect.top + y};
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('click', o));
  }, {vn: viewName, col: column, x: CAT_NAME_X,
    y: CAT_ROW_TOP + CAT_ROW_PITCH * rowIndex + CAT_ROW_CENTRE});
}

async function openAgeCardMenu(page: import('@playwright/test').Page): Promise<void> {
  const cards = ageCard(page);
  expect(await cards.count(),
    'the laid-out Filter Panel does not hold exactly one card captioned AGE — the menu below would be opened on an ambiguous card')
    .toBe(1);
  const card = cards.first();
  const cardBox = await card.boundingBox();
  if (!cardBox) throw new Error('the AGE card in the laid-out Filter Panel has no layout box');
  await page.mouse.move(cardBox.x + cardBox.width / 2, cardBox.y + cardBox.height / 2, {steps: 6});
  const indicator = card.locator('.d4-filter-indicator').first();
  await expect(indicator,
    'the AGE card\'s filter-options indicator did not become visible under the pointer — it carries visibility:hidden until the card is hovered, so the menu cannot be opened')
    .toBeVisible({timeout: 10_000});
  const iBox = await indicator.boundingBox();
  if (!iBox) throw new Error('the AGE card\'s filter-options indicator reports no layout box to click');
  await page.mouse.click(iBox.x + iBox.width / 2, iBox.y + iBox.height / 2);
  await expect.poll(async () => page.locator('.d4-menu-popup').count(), {
    message: 'clicking the AGE card\'s indicator did not leave exactly one .d4-menu-popup open',
    timeout: 10_000, intervals: [200, 400, 800],
  }).toBe(1);
  await expect(page.locator('.d4-menu-popup [name="div-Missing-values"]').first(),
    'the open card menu carries no [name="div-Missing-values"] group')
    .toBeAttached({timeout: 10_000});
}

async function missingValuesGlyphs(
  page: import('@playwright/test').Page,
): Promise<Record<string, string>> {
  return page.evaluate((leaves: [string, string][]) => {
    const popups = document.querySelectorAll('.d4-menu-popup');
    if (popups.length !== 1)
      throw new Error(`expected exactly one .d4-menu-popup while reading the radio glyphs, found ${popups.length}`);
    const out: Record<string, string> = {};
    for (const [key, name] of leaves) {
      const leaf = popups[0].querySelector(`[name="${name}"]`);
      if (!leaf)
        throw new Error(`the Missing values leaf [name="${name}"] is not present in the open card menu`);
      const glyph = leaf.querySelector('.d4-menu-item-check > i');
      if (!glyph)
        throw new Error(`the Missing values leaf [name="${name}"] carries no .d4-menu-item-check > i radio glyph`);
      out[key] = glyph.getAttribute('name') ?? '';
    }
    return out;
  }, MV_LEAVES);
}

async function missingValuesGroupDisplay(page: import('@playwright/test').Page): Promise<string> {
  return page.evaluate(() => {
    const group = document.querySelector('.d4-menu-popup [name="div-Missing-values"]');
    if (!group) throw new Error('no [name="div-Missing-values"] group in the open card menu');
    const container = group.querySelector('.d4-menu-item-container');
    if (!container)
      throw new Error('the Missing values group carries no nested .d4-menu-item-container');
    return getComputedStyle(container).display;
  });
}

async function inertChromePoint(
  page: import('@playwright/test').Page,
): Promise<{x: number; y: number}> {
  const pt = await page.evaluate(() => {
    const w = window.innerWidth;
    for (let y = 6; y < 78; y += 6) {
      for (let x = w - 8; x > w * 0.35; x -= 14) {
        const el = document.elementFromPoint(x, y) as HTMLElement | null;
        if (!el || el.tagName !== 'DIV' || el.hasAttribute('name') || el.children.length === 0) continue;
        if (el.closest('.d4-menu-popup') || el.closest('.d4-filter') || el.closest('.d4-dialog')) continue;
        if (el.querySelector('canvas') || el.querySelector('input')) continue;
        const r = el.getBoundingClientRect();
        if (r.width < 100) continue;
        return {x, y};
      }
    }
    return null;
  });
  if (!pt) {
    throw new Error('no inert chrome point was found to dismiss the card menu with; Escape must not be '
      + 'substituted here — it resets the card, unchecks it and pushes the row count back to the full table');
  }
  return pt;
}

async function dismissCardMenu(page: import('@playwright/test').Page): Promise<void> {
  const pt = await inertChromePoint(page);
  await page.mouse.click(pt.x, pt.y);
  await expect.poll(async () => page.locator('.d4-menu-popup').count(), {
    message: `clicking inert chrome at (${pt.x}, ${pt.y}) did not close the card menu`,
    timeout: 10_000, intervals: [200, 400, 800],
  }).toBe(0);
}

async function laidOutPanelCaptions(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() => {
    const panels = Array.from(document.querySelectorAll('[name="viewer-Filters"]'))
      .filter((e) => (e as HTMLElement).offsetParent !== null);
    if (panels.length !== 1)
      throw new Error(`expected exactly one laid-out Filter Panel to read, found ${panels.length}`);
    return Array.from(panels[0].querySelectorAll('.d4-filter .d4-filter-column-name'))
      .map((e) => (e.textContent ?? '').trim());
  });
}

interface ComboPick { opened: boolean; accepted: string; added: string[]; }

async function drivePanelColumnCombo(
  page: import('@playwright/test').Page, column: string,
): Promise<ComboPick> {
  const before = await laidOutPanelCaptions(page);
  await page.evaluate(() => {
    const panels = Array.from(document.querySelectorAll('[name="viewer-Filters"]'))
      .filter((e) => (e as HTMLElement).offsetParent !== null);
    if (panels.length !== 1)
      throw new Error(`expected exactly one laid-out Filter Panel to drive, found ${panels.length}`);
    const combo = panels[0].querySelector('.d4-filter-group-header [name="div-column-combobox-"]');
    if (!combo)
      throw new Error('the laid-out Filter Panel header exposes no [name="div-column-combobox-"]');
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const label = combo.querySelector('.d4-column-selector-column');
    (label ?? combo).dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  });
  const opened = await page.waitForFunction(
    () => !!document.querySelector('.d4-column-selector-backdrop'), null, {timeout: 8000})
    .then(() => true).catch(() => false);
  if (!opened) return {opened: false, accepted: '', added: []};
  await page.keyboard.press(column[0].toLowerCase());
  const search = page.locator('input.d4-column-selector-search-input').first();
  await expect(search,
    `the open column picker never produced its search box after "${column[0]}" was typed at it, so the pick `
    + 'could not be addressed by name and Enter would commit whichever row the picker happened to highlight')
    .toBeAttached({timeout: 10_000});
  await page.keyboard.press('Control+a');
  await page.keyboard.type(column, {delay: 40});
  const accepted = await search.inputValue();
  await page.keyboard.press('Enter');
  await expect.poll(async () => page.locator('.d4-column-selector-backdrop').count(), {
    message: `the column picker never closed after "${column}" was committed with Enter — the pick did not land`,
    timeout: 15_000, intervals: [300, 600, 1200],
  }).toBe(0);
  const after = await laidOutPanelCaptions(page);
  const rest = [...before];
  const added: string[] = [];
  for (const caption of after) {
    const at = rest.indexOf(caption);
    if (at === -1) added.push(caption);
    else rest.splice(at, 1);
  }
  return {opened: true, accepted, added};
}

test('Filters — Cloned View Synchronization', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});
  await installViewResolver(page);

  const ORIG = 'demog';
  let CLONE = '';

  await page.evaluate((name: string) => {
    const tv = grok.shell.tv;
    tv.dataFrame.name = name;
    tv.name = name;
  }, ORIG);
  await v.resetFilters(page);
  const initial = await page.evaluate((name: string) => {
    const tv = (window as any).__tv(name);
    return {
      viewName: tv.name,
      rowCount: tv.dataFrame.rowCount,
      trueCount: tv.dataFrame.filter.trueCount,
      cards: ((window as any).__cards(name) as HTMLElement[]).length,
    };
  }, ORIG);
  expect(initial.viewName).toBe(ORIG);
  expect(initial.rowCount).toBe(FULL);
  expect(initial.trueCount).toBe(FULL);
  expect(initial.cards).toBe(0);

  const cats: {race: string[]; sex: string[]} = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {race: [...df.col('RACE').categories], sex: [...df.col('SEX').categories]};
  });
  expect(cats.race.filter((c) => c !== 'Asian').length,
    'the RACE column holds no category other than "Asian", so the setup criterion would select nothing')
    .toBeGreaterThan(0);
  await v.applyCategoricalFilter(page, 'RACE', cats.race.filter((c) => c !== 'Asian'));
  await v.applyCategoricalFilter(page, 'SEX', cats.sex);
  const sexFirst = cats.sex.filter((c) => c != null)[0];
  expect(sexFirst, 'the SEX column offers no non-empty category to narrow to').toBeTruthy();
  const truePanel = await trueCountOf(page, ORIG);
  expect(truePanel).toBeLessThan(FULL);
  expect(truePanel).toBeGreaterThan(0);
  const setupFilters = await readFilters(page, ORIG);
  expect(setupFilters.map((f) => f.col).sort()).toEqual(['RACE', 'SEX']);
  expect(setupFilters.filter((f) => f.isFiltering).map((f) => f.col)).toEqual(['RACE']);

  try {
    await softStep('Scenario 1 - Step 5 clone inherits card set and per-card isActive/isFiltering states', async () => {
      const origCaptions = await panelCaptionsOf(page, ORIG);
      expect(origCaptions.length).toBeGreaterThan(0);
      const origFilters = await readFilters(page, ORIG);
      const namesBefore = await viewNames(page);

      expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Clone View']),
        'the View | Layout | Clone View menu leaf was not actuated').toBe(true);
      await expect.poll(async () => (await viewNames(page)).length,
        {timeout: 15_000, intervals: [400, 800, 1500]}).toBe(namesBefore.length + 1);
      const namesAfter = await viewNames(page);
      const addedNames = namesAfter.filter((n) => !namesBefore.includes(n));
      expect(addedNames,
        `Clone View must add exactly one distinctly-named TableView (before: ${namesBefore.join(' | ')}; after: ${namesAfter.join(' | ')})`)
        .toHaveLength(1);
      CLONE = addedNames[0];

      await expect.poll(async () => page.evaluate((vn: string) =>
        ((window as any).__tv(vn).root as HTMLElement).querySelectorAll('[name="viewer-Filters"]').length, CLONE),
      {timeout: 10_000, intervals: [400, 800, 1500]}).toBeGreaterThan(0);
      expect(await page.evaluate((vn: string) => {
        const el = ((window as any).__tv(vn).root as HTMLElement)
          .querySelector('[name="viewer-Filters"]') as HTMLElement;
        return el != null && el.offsetParent !== null;
      }, CLONE)).toBe(true);

      const cloneCaptions = await panelCaptionsOf(page, CLONE);
      expect(cloneCaptions).toEqual(origCaptions);

      const cloneFilters = await readFilters(page, CLONE);
      const key = (f: {col: string | null; type: string; isActive: boolean; isFiltering: boolean}) =>
        `${f.col}|${f.type}|${f.isActive}|${f.isFiltering}`;
      const origKeys = origFilters.filter((f) => f.col).map(key).sort();
      const cloneKeys = cloneFilters.filter((f) => f.col).map(key).sort();
      expect(origKeys.length,
        'not one filter in the ORIGINAL reported a column name — every read fell into the proxy catch '
        + 'and returned col: null, so the two-sided comparison below would be [] vs [] and could not fail')
        .toBeGreaterThan(0);
      expect(cloneKeys.length,
        'not one filter in the CLONE reported a column name — the comparison below would be satisfied '
        + 'by a clone whose filter group is empty')
        .toBe(origKeys.length);
      expect(cloneKeys).toEqual(origKeys);

      expect(cloneFilters.filter((f) => f.isFiltering).map((f) => f.col)).toEqual(['RACE']);
    });

    await softStep('Scenario 1 - Step 7 change RACE in original through its card, clone criterion mirrors', async () => {
      const before = await trueCountOf(page, ORIG);
      const origRaceBefore = await selectedOf(page, ORIG, 'RACE');
      const cloneRaceBefore = await selectedOf(page, CLONE, 'RACE');
      expect(cloneRaceBefore.length,
        `the clone's RACE criterion is already a single category (${JSON.stringify(cloneRaceBefore)}) before `
        + 'the change, so reading a single category back from it afterwards would be a pre-satisfied assertion')
        .toBeGreaterThan(1);
      expect(cloneRaceBefore,
        `the clone's RACE criterion (${JSON.stringify(cloneRaceBefore)}) is not the original's `
        + `(${JSON.stringify(origRaceBefore)}) before the change, so what the clone reads afterwards would not `
        + 'be a move away from a shared starting point')
        .toEqual(origRaceBefore);

      await clickCategoryNameRowIn(page, ORIG, 'RACE', 2);
      await expect.poll(async () => (await selectedOf(page, ORIG, 'RACE')).length,
        {message: 'clicking a category name row in the RACE card of the ORIGINAL view did not narrow that card to a '
          + 'single category — the exclusive-select gesture never landed, so nothing was changed for the '
          + 'clone to follow',
        timeout: 20_000, intervals: [300, 600, 1200]}).toBe(1);
      const origRaceAfter = await selectedOf(page, ORIG, 'RACE');
      const narrowed = await trueCountOf(page, ORIG);
      expect(narrowed).not.toBe(before);
      expect(narrowed).toBeGreaterThan(0);

      await expect.poll(async () => JSON.stringify(await selectedOf(page, CLONE, 'RACE')),
        {message: 'the RACE card of the CLONE did not take on the criterion the ORIGINAL card was clicked '
          + `to (${JSON.stringify(origRaceAfter)}); it started from the wider `
          + `${JSON.stringify(cloneRaceBefore)}. A criterion set through a card gesture is broadcast on the `
          + 'DataFrame event bus and copied into every other card on the same (column, filter type) — '
          + 'grid_filter_base.dart:111-125 fires FILTER_CRITERIA_CHANGED, grid_filter_base.dart:71-81 '
          + 'applies the saveState() of the sender onto this card — so the clone had to reach that same set',
        timeout: 20_000, intervals: [300, 600, 1200]}).toBe(JSON.stringify(origRaceAfter));
      const cloneFiltering = (await readFilters(page, CLONE)).filter((f) => f.isFiltering);
      expect(cloneFiltering.map((f) => f.col)).toEqual(['RACE']);
    });

    await softStep('Scenario 1 - Step 9 disable SEX in original mirrors to clone (isActive mirror)', async () => {
      const beforeSex = await page.evaluate(({o, c}) => {
        const find = (n: string) => (window as any).__tv(n).getFiltersGroup().filters
          .find((f: any) => { try { return f.filterColumnName === 'SEX' && f.filterType === 'categorical'; } catch (_) { return false; } });
        return {orig: find(o)?.isActive, clone: find(c)?.isActive};
      }, {o: ORIG, c: CLONE});
      expect(beforeSex).toEqual({orig: true, clone: true});

      const tcBeforeNarrow = await trueCountOf(page, ORIG);
      await activateView(page, ORIG);
      const narrowedCount = (await v.applyCategoricalFilter(page, 'SEX', [sexFirst])).filteredCount;
      const narrowedFiltering = await page.evaluate((vn: string) => {
        const sex = (window as any).__tv(vn).getFiltersGroup().filters
          .find((f: any) => { try { return f.filterColumnName === 'SEX' && f.filterType === 'categorical'; } catch (_) { return false; } });
        if (!sex) throw new Error(`view "${vn}" has no SEX categorical card to narrow`);
        return sex.isFiltering;
      }, ORIG);
      expect(narrowedFiltering,
        'the SEX card is not filtering after being narrowed to one category — the disable below would have nothing to lift')
        .toBe(true);
      expect(narrowedCount,
        'narrowing SEX to one category did not drop any rows').toBeLessThan(tcBeforeNarrow);
      expect(narrowedCount, 'narrowing SEX to one category emptied the table').toBeGreaterThan(0);

      await clickCardCheckboxIn(page, ORIG, 'SEX');
      const isActiveOfSex = async () => page.evaluate(({o, c}) => {
        const sexOf = (n: string) => (window as any).__tv(n).getFiltersGroup().filters
          .find((f: any) => { try { return f.filterColumnName === 'SEX' && f.filterType === 'categorical'; } catch (_) { return false; } });
        return {orig: sexOf(o)?.isActive, clone: sexOf(c)?.isActive};
      }, {o: ORIG, c: CLONE});
      await expect.poll(isActiveOfSex,
        {message: 'unticking the SEX card\'s own checkbox in the ORIGINAL did not switch that card off in '
          + 'both views — the isActive mirror never crossed to the clone\'s filter group',
        timeout: 20_000, intervals: [400, 800, 1500]}).toEqual({orig: false, clone: false});
      const tcAfterDisable = await trueCountOf(page, ORIG);
      expect(tcAfterDisable,
        'disabling the SEX card did not restore the rows its criterion had removed').toBe(tcBeforeNarrow);
      expect(tcAfterDisable,
        'disabling the SEX card left the row set where the narrowing had put it').toBeGreaterThan(narrowedCount);
      expect(await cardDisabledIn(page, ORIG, 'SEX'),
        'the SEX card in the ORIGINAL view is not painted disabled').toBe(true);
      expect(await cardDisabledIn(page, CLONE, 'SEX'),
        'the SEX card in the CLONE view is not painted disabled').toBe(true);
      expect(await cardDisabledIn(page, ORIG, 'RACE')).toBe(false);
      expect(await cardDisabledIn(page, CLONE, 'RACE')).toBe(false);
    });

    await softStep('Scenario 1 - Step 10 different-type card on same column is not mirrored', async () => {
      const tcBeforeHist = await trueCountOf(page, ORIG);

      await activateView(page, ORIG);
      const tcWithHist = await v.applyNumericFilter(page, 'AGE', 30, 60);
      const hist = await page.evaluate((vn: string) => {
        const st = (window as any).__tv(vn).getFiltersGroup().getStates('AGE', 'histogram');
        if (!st || !st.length)
          throw new Error(`AGE histogram card was not created in view "${vn}"`);
        return {count: st.length, active: st[0].active};
      }, ORIG);
      expect(hist.count).toBe(1);
      expect(hist.active).toBe(true);
      expect(tcWithHist).toBeLessThan(tcBeforeHist);
      expect(tcWithHist).toBeGreaterThan(0);

      const added = await page.evaluate(async (vn: string) => {
        const fg = (window as any).__tv(vn).getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'AGE'});
        await new Promise((r) => setTimeout(r, 1500));
        const cat = fg.filters.find((f: any) => { try { return f.filterColumnName === 'AGE' && f.filterType === 'categorical'; } catch (_) { return false; } });
        if (!cat)
          throw new Error(`AGE categorical card was not created in view "${vn}"`);
        return {active: cat.isActive};
      }, CLONE);
      expect(added.active).toBe(true);
      const tcBeforeToggle = await trueCountOf(page, ORIG);
      expect(tcBeforeToggle).toBeGreaterThan(0);

      const cloneAgeIsActive = async () => page.evaluate((vn: string) => {
        const cat = (window as any).__tv(vn).getFiltersGroup().filters
          .find((f: any) => { try { return f.filterColumnName === 'AGE' && f.filterType === 'categorical'; } catch (_) { return false; } });
        if (!cat)
          throw new Error(`AGE categorical card vanished from view "${vn}" before the toggle`);
        return cat.isActive;
      }, CLONE);
      expect(await cloneAgeIsActive(),
        'the clone\'s AGE card is already switched off, so switching it off proves nothing').toBe(true);
      await clickCardCheckboxIn(page, CLONE, 'AGE');
      await expect.poll(cloneAgeIsActive,
        {message: 'unticking the AGE card\'s own checkbox in the CLONE did not switch that card off',
          timeout: 20_000, intervals: [400, 800, 1500]}).toBe(false);

      const origAfter = await page.evaluate((vn: string) => {
        const st = (window as any).__tv(vn).getFiltersGroup().getStates('AGE', 'histogram');
        if (!st || !st.length)
          throw new Error(`AGE histogram card disappeared from view "${vn}"`);
        return {count: st.length, active: st[0].active};
      }, ORIG);
      expect(origAfter.count).toBe(1);
      expect(origAfter.active).toBe(true);

      const origSamples: string[] = [];
      for (let i = 0; i < 6; i++) {
        await page.waitForTimeout(400);
        origSamples.push(JSON.stringify(await page.evaluate((vn: string) => {
          const view = (window as any).__tv(vn);
          const st = view.getFiltersGroup().getStates('AGE', 'histogram');
          return {count: st ? st.length : 0, active: st && st.length ? st[0].active : null,
            trueCount: view.dataFrame.filter.trueCount};
        }, ORIG)));
      }
      expect(Array.from(new Set(origSamples)),
        'the ORIGINAL\'s AGE histogram card and the narrowing it holds had to stay put at every sample of '
        + 'a 2.4s window after the clone-side toggle — a mirror that crosses filter types on a longer '
        + `debounce would otherwise land after a single snapshot was taken; samples: ${origSamples.join(' ; ')}`)
        .toEqual([JSON.stringify({count: 1, active: true, trueCount: tcBeforeToggle})]);
      const tcAfter = await trueCountOf(page, ORIG);
      expect(tcAfter).toBe(tcBeforeToggle);
      expect(tcAfter).toBeLessThan(tcBeforeHist);
    });

    await softStep('Scenario 1 - Step 11 removing the clone\'s card leaves the original card and its narrowing', async () => {
      await clickCardCheckboxIn(page, ORIG, 'SEX');
      await expect.poll(async () => (await readFilters(page, ORIG))
        .find((f) => f.col === 'SEX' && f.type === 'categorical')?.isActive,
      {message: 'ticking the SEX card\'s own checkbox back on in the original did not re-enable it',
        timeout: 20_000, intervals: [400, 800, 1500]}).toBe(true);
      await activateView(page, ORIG);
      await v.applyCategoricalFilter(page, 'SEX', [sexFirst]);
      const armedCategory = String(sexFirst);
      await expect.poll(async () => (await readFilters(page, ORIG))
        .find((f) => f.col === 'SEX' && f.type === 'categorical')?.isFiltering,
      {message: 'the SEX card in the original must be filtering before the removal',
        timeout: 20_000, intervals: [400, 800, 1500]}).toBe(true);
      const tcBeforeRemoval = await trueCountOf(page, ORIG);

      await activateView(page, CLONE);
      await page.evaluate(({vn, col}) => {
        const x = (window as any).__card(vn, col).querySelector('[name="icon-times"]') as HTMLElement | null;
        if (!x)
          throw new Error(`the "${col}" card in view "${vn}" carries no [name="icon-times"] remove icon — the removal gesture never happened`);
        x.click();
      }, {vn: CLONE, col: 'SEX'});
      await expect.poll(async () => page.evaluate((vn: string) =>
        ((window as any).__cards(vn) as HTMLElement[])
          .map((e) => e.querySelector('.d4-filter-column-name')?.textContent?.trim()), CLONE),
      {message: 'the SEX card is still in the clone — the X click did not remove it',
        timeout: 20_000, intervals: [400, 800, 1500]}).not.toContain('SEX');

      expect(await panelCaptionsOf(page, ORIG),
        'removing the card in the clone also removed it from the original').toContain('SEX');
      const origSex = (await readFilters(page, ORIG))
        .find((f) => f.col === 'SEX' && f.type === 'categorical');
      expect(origSex, 'the original has no SEX categorical state left after the clone-side removal').toBeTruthy();
      expect(origSex!.isActive, 'the original\'s SEX card was switched off by the clone-side removal').toBe(true);
      expect(origSex!.isFiltering, 'the original\'s SEX card stopped filtering after the clone-side removal').toBe(true);
      expect(await selectedOf(page, ORIG, 'SEX'),
        'the original\'s SEX criterion was rewritten by the clone-side removal').toEqual([armedCategory]);
      expect(await trueCountOf(page, ORIG),
        'the filtered row count moved when a card was removed in the other view').toBe(tcBeforeRemoval);
    });

    await softStep('Scenario 2 - Step 4 missing-values choice set from the AGE card menu survives the clone', async () => {
      await v.openTable(page, {path: datasetPath, withFilterPanel: true});
      const MV: string = await page.evaluate(() => grok.shell.tv.name);
      const base = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return {rowCount: df.rowCount, missing: df.col('AGE').stats.missingValueCount,
          trueCount: df.filter.trueCount};
      });
      expect(base.rowCount,
        `the missing-values check needs the 5850-row demog golden table; got ${base.rowCount}`).toBe(FULL);
      expect(base.missing,
        `AGE must hold exactly one missing value for the ${FULL} -> ${FULL - 1} delta below to be the missing-value row; got ${base.missing}`)
        .toBe(1);
      expect(base.trueCount,
        `nothing may be filtered before the choice is made, otherwise the delta is not the option's; trueCount ${base.trueCount}`)
        .toBe(FULL);
      const origCounts = await panelCounts(page, MV);
      expect(origCounts.ageCards,
        `the default Filter Panel must paint exactly one AGE card to drive; got ${origCounts.ageCards} of ${origCounts.cards} cards`)
        .toBe(1);

      await openAgeCardMenu(page);
      expect(await missingValuesGlyphs(page),
        'the AGE card\'s Missing values group does not start on "Keep missing values" with the other two clear — without that baseline the read-back after the clone proves nothing')
        .toEqual({keep: DOT, filterOut: CIRCLE, showOnly: CIRCLE});

      expect(await missingValuesGroupDisplay(page),
        'the Missing values group is already expanded before it is hovered — the hover below could not be shown to be what opened it')
        .toBe('none');
      const groupBox = await boxOf(page, '.d4-menu-popup [name="div-Missing-values"]');
      await page.mouse.move(groupBox.x + groupBox.width / 2, groupBox.y + groupBox.height / 2, {steps: 6});
      await expect.poll(async () => missingValuesGroupDisplay(page), {
        message: 'hovering the Missing values group did not expand it — its nested .d4-menu-item-container never left display:none',
        timeout: 10_000, intervals: [200, 400, 800],
      }).toBe('flex');

      const leafBox = await boxOf(page, '.d4-menu-popup [name="div-Missing-values---Filter-out-missing-values"]');
      await page.mouse.click(leafBox.x + leafBox.width / 2, leafBox.y + leafBox.height / 2);
      await expect.poll(async () => trueCountOf(page, MV), {
        message: `choosing "Filter out missing values" on the AGE card did not remove the single missing-value row (expected ${FULL - 1})`,
        timeout: 20_000, intervals: [300, 600, 1200],
      }).toBe(FULL - 1);

      await dismissCardMenu(page);
      expect(await trueCountOf(page, MV),
        `dismissing the card menu moved the filtered row count — the dismissal was not inert (expected ${FULL - 1})`)
        .toBe(FULL - 1);

      await openAgeCardMenu(page);
      expect(await missingValuesGlyphs(page),
        'the AGE card menu does not read "Filter out missing values" back as the selected option immediately after it was chosen there')
        .toEqual({keep: CIRCLE, filterOut: DOT, showOnly: CIRCLE});
      await dismissCardMenu(page);

      const mvNamesBefore = await viewNames(page);
      expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Clone View']),
        'the View | Layout | Clone View menu leaf was not actuated').toBe(true);
      await expect.poll(async () => (await viewNames(page)).length, {
        message: 'Clone View did not add a second TableView',
        timeout: 15_000, intervals: [400, 800, 1500],
      }).toBe(mvNamesBefore.length + 1);
      const mvNamesAfter = await viewNames(page);
      const mvAdded = mvNamesAfter.filter((n) => !mvNamesBefore.includes(n));
      expect(mvAdded,
        `Clone View must add exactly one distinctly-named TableView (before: ${mvNamesBefore.join(' | ')}; after: ${mvNamesAfter.join(' | ')})`)
        .toHaveLength(1);
      const MV_CLONE = mvAdded[0];
      expect(await trueCountOf(page, MV),
        `cloning the view moved the filtered row count (expected ${FULL - 1})`).toBe(FULL - 1);

      await activateView(page, MV_CLONE);
      const cloneBefore = await panelCounts(page, MV_CLONE);
      expect(cloneBefore.ageCards,
        `the clone must come up with exactly one AGE card before the re-add is attempted; got ${cloneBefore.ageCards} of ${cloneBefore.cards} cards`)
        .toBe(1);
      const cardless: string = await page.evaluate((vn: string) => {
        const painted = new Set(((window as any).__cards(vn) as HTMLElement[])
          .map((e) => e.querySelector('.d4-filter-column-name')?.textContent?.trim()));
        const free = (window as any).__tv(vn).dataFrame.columns.toList()
          .filter((c: any) => (c.isCategorical || c.isNumerical) && !painted.has(c.name))
          .map((c: any) => c.name);
        return free[free.length - 1] ?? '';
      }, MV_CLONE);
      expect(cardless,
        'every column of the table already has a card in the clone\'s panel, so the picker has nothing it '
        + 'could legitimately add and the positive control below cannot be run')
        .toBeTruthy();
      const positive = await drivePanelColumnCombo(page, cardless);
      expect(positive.opened,
        `the panel column picker never opened in the clone for the card-less column "${cardless}"`)
        .toBe(true);
      expect(positive.accepted,
        `the clone's column picker did not hold "${cardless}" when Enter committed the pick — the picker commits `
        + 'the name its search box carries (column_combo_box.dart:361-364), so a commit made against any other '
        + `text is a pick of some other column (the picker read "${positive.accepted}")`)
        .toBe(cardless);
      await expect.poll(async () => (await panelCounts(page, MV_CLONE)).cards,
        {message: `picking the card-less column "${cardless}" in the clone's column selector added no card — `
          + 'the picker adds nothing for ANY column, so the refusal asserted below for AGE would be the '
          + 'picker being broken rather than the product declining a second instance on the same column',
        timeout: 20_000, intervals: [400, 800, 1500]}).toBe(cloneBefore.cards + 1);
      expect(positive.added,
        `the card the picker added is not captioned "${cardless}" — the column is deliberately the LAST card-free `
        + 'one rather than the first, so a picker that ignored the typed name and committed whichever row it had '
        + `highlighted would add some other column's card here and the AGE refusal below would be measured against `
        + `a picker that never accepts what it is told (added: ${positive.added.join(', ') || '(nothing)'})`)
        .toEqual([cardless]);
      const cloneWithExtra = await panelCounts(page, MV_CLONE);
      expect(cloneWithExtra.ageCards,
        `adding a card for "${cardless}" also touched the AGE cards`).toBe(1);
      expect(await trueCountOf(page, MV),
        `adding the "${cardless}" card moved the filtered row count, so the reads that follow are no longer about the missing-values choice`)
        .toBe(FULL - 1);

      const refusal = await drivePanelColumnCombo(page, 'AGE');
      expect(refusal.opened,
        'the panel column picker never opened in the clone, so the product\'s refusal to add a second AGE card was never exercised')
        .toBe(true);
      expect(refusal.accepted,
        'the clone\'s column picker did not hold "AGE" when Enter committed the pick, so AGE was never the column '
        + 'the picker was asked for and the card set holding still below says nothing about a second AGE instance '
        + `(the picker read "${refusal.accepted}")`)
        .toBe('AGE');
      expect(refusal.added,
        `picking AGE added a card captioned ${refusal.added.join(', ')} — the pick landed on some other column`)
        .toEqual([]);
      const cloneSamples: string[] = [];
      for (let i = 0; i < 8; i++) {
        await page.waitForTimeout(500);
        cloneSamples.push(JSON.stringify(await panelCounts(page, MV_CLONE)));
      }
      expect(Array.from(new Set(cloneSamples)),
        'picking AGE in the clone\'s column selector changed the clone\'s card set over a sustained 4s window; '
        + 'addDefaultFilter (filters_core.dart:829-841) dedupes on (column, default filter type) and only scrolls '
        + `the existing card into view, so every sample must read ${JSON.stringify(cloneWithExtra)}; samples: ${cloneSamples.join(' ; ')}`)
        .toEqual([JSON.stringify(cloneWithExtra)]);

      await activateView(page, MV);
      expect(await trueCountOf(page, MV),
        `switching back to the original view lost the missing-values narrowing (expected ${FULL - 1})`)
        .toBe(FULL - 1);
      await openAgeCardMenu(page);
      const finalGlyphs = await missingValuesGlyphs(page);
      const finalTrueCount = await trueCountOf(page, MV);
      expect(finalGlyphs,
        `github-1984: the original view's AGE card no longer reads "Filter out missing values" in its own menu after the clone round trip (glyphs ${JSON.stringify(finalGlyphs)}, trueCount ${finalTrueCount})`)
        .toEqual({keep: CIRCLE, filterOut: DOT, showOnly: CIRCLE});
      expect(finalTrueCount,
        `github-1984: the original view's filtered row count moved during the clone round trip (expected ${FULL - 1}, glyphs ${JSON.stringify(finalGlyphs)})`)
        .toBe(FULL - 1);
      await dismissCardMenu(page);
    });

    await softStep('Scenario 3 - Step 4 filtering the host frame leaves the shared-column frame untouched', async () => {
      const before = await page.evaluate(async () => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 800));
        const t1 = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
        const v1 = grok.shell.addTableView(t1);
        v1.name = 'HostView';
        await new Promise((r) => setTimeout(r, 1200));
        const sexCol = t1.col('SEX');
        const t2 = DG.DataFrame.fromColumns([sexCol]);
        grok.shell.addTableView(t2).name = 'SharedColView';
        await new Promise((r) => setTimeout(r, 1200));

        t1.col('SEX').setTag('sharedColumnProbe', 'shared');
        const sharesColumnObject = t2.col('SEX').getTag('sharedColumnProbe') === 'shared';

        grok.shell.v = v1;
        v1.getFiltersGroup();
        return {t1: t1.filter.trueCount, t2: t2.filter.trueCount, t2rows: t2.rowCount,
          sharesColumnObject, firstCat: String(t1.col('SEX').categories[0])};
      });
      await activateView(page, 'HostView');
      await v.applyCategoricalFilter(page, 'SEX', [before.firstCat]);
      const samples: {t1: number; t2: number}[] = [];
      for (let i = 0; i < 8; i++) {
        if (i > 0) await page.waitForTimeout(500);
        samples.push(await page.evaluate(() => ({
          t1: (window as any).__tv('HostView').dataFrame.filter.trueCount,
          t2: (window as any).__tv('SharedColView').dataFrame.filter.trueCount,
        })));
      }
      const outcome = {before, after: samples[samples.length - 1]};
      expect(outcome.before.sharesColumnObject,
        't2 does not share the SEX column object with t1 — the isolation check below could not fail')
        .toBe(true);
      expect(outcome.before.t1).toBe(FULL);
      expect(outcome.before.t2).toBe(outcome.before.t2rows);
      expect(outcome.after.t1).toBeLessThan(FULL);
      expect(outcome.after.t1).toBeGreaterThan(0);
      expect(Array.from(new Set(samples.map((s) => s.t2))),
        'the frame that merely shares the SEX column object moved at some point in a sustained 3.5s window after '
        + `the host frame was filtered, so its isolation is not a snapshot taken before contamination arrived `
        + `(every sample had to read ${outcome.before.t2rows}; samples: ${samples.map((s) => s.t2).join(', ')})`)
        .toEqual([outcome.before.t2rows]);
    });

    await softStep('Scenario 4 - Step 3 saved layout restores panel, cards and trueCount', async () => {
      let layoutId: string | null = null;
      const LAYOUT_VIEW = 'LayoutHost';
      try {
        await page.evaluate(async (vn: string) => {
          grok.shell.closeAll();
          await new Promise((r) => setTimeout(r, 800));
          const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
          const view = grok.shell.addTableView(df);
          df.name = vn;
          view.name = vn;
          await new Promise((r) => setTimeout(r, 1500));
          view.getFiltersGroup();
        }, LAYOUT_VIEW);
        await v.resetFilters(page);
        await activateView(page, LAYOUT_VIEW);
        await v.applyCategoricalFilter(page, 'RACE', ['Asian']);
        const withBothCards = await v.applyNumericFilter(page, 'AGE', 30, 60);

        await clickCardCheckboxIn(page, LAYOUT_VIEW, 'RACE');
        await expect.poll(async () => page.evaluate((vn: string) =>
          (window as any).__tv(vn).getFiltersGroup().getStates('RACE', 'categorical')[0]?.active, LAYOUT_VIEW),
        {message: 'unticking the RACE card\'s own checkbox did not switch that card off, so the layout '
          + 'about to be saved does not carry the one-card-off state the restore is measured against',
        timeout: 20_000, intervals: [400, 800, 1500]}).toBe(false);

        const saved = await page.evaluate(async (vn: string) => {
          const view = (window as any).__tv(vn);
          const fg = view.getFiltersGroup();
          const raceEl = (window as any).__card(vn, 'RACE') as HTMLElement;
          const ageEl = (window as any).__card(vn, 'AGE') as HTMLElement;
          const layout = view.saveLayout();
          await grok.dapi.layouts.save(layout);
          let serverFound = false;
          for (let i = 0; i < 20 && !serverFound; i++) {
            try { serverFound = !!(await grok.dapi.layouts.find(layout.id)); }
            catch (_) { serverFound = false; }
            if (!serverFound) await new Promise((r) => setTimeout(r, 500));
          }
          return {
            serverFound,
            savedTrueCount: view.dataFrame.filter.trueCount,
            captions: ((window as any).__cards(vn) as HTMLElement[])
              .map((e) => (e.querySelector('.d4-filter-column-name')?.textContent ?? '').trim()),
            layoutId: layout.id,
            raceActive: fg.getStates('RACE', 'categorical')[0]?.active,
            ageActive: fg.getStates('AGE', 'histogram')[0]?.active,
            raceDisabledClass: raceEl.classList.contains('d4-filter-disabled'),
            ageDisabledClass: ageEl.classList.contains('d4-filter-disabled'),
          };
        }, LAYOUT_VIEW);
        layoutId = saved.layoutId;
        expect(saved.serverFound,
          'the saved layout cannot be fetched back from the server — saveLayout() stamps the id '
          + 'client-side before the round-trip, so the re-apply below would restore nothing')
          .toBe(true);
        expect(saved.raceActive).toBe(false);
        expect(saved.raceDisabledClass).toBe(true);
        expect(saved.ageActive).toBe(true);
        expect(saved.ageDisabledClass).toBe(false);
        expect(saved.savedTrueCount).toBeGreaterThan(withBothCards);
        expect(saved.savedTrueCount).toBeGreaterThan(0);
        expect(saved.savedTrueCount).toBeLessThan(FULL);
        expect(saved.captions.length).toBeGreaterThan(0);

        await v.applyNumericFilter(page, 'AGE', 0, 200);
        await page.evaluate(async (vn: string) => {
          const view = (window as any).__tv(vn);
          view.getFiltersGroup().close();
          await new Promise((r) => setTimeout(r, 800));
          view.addViewer('Bar chart');
          await new Promise((r) => setTimeout(r, 1500));
        }, LAYOUT_VIEW);
        expect(await page.locator('[name="viewer-Filters"]').count()).toBe(0);
        const perturbed = await trueCountOf(page, LAYOUT_VIEW);
        expect(perturbed,
          'the perturbation left df.filter.trueCount at the saved value — the restore assert below could not fail')
          .not.toBe(saved.savedTrueCount);

        const restored = await page.evaluate(async ({id, vn, cap}) => {
          const view = (window as any).__tv(vn);
          const s = await grok.dapi.layouts.find(id);
          if (!s) throw new Error(`the saved layout ${id} is not on the server, so nothing can be re-applied`);
          const applied = new Promise<void>((resolve) => {
            const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
            setTimeout(resolve, cap);
          });
          view.loadLayout(s);
          await applied;
          const root = view.root as HTMLElement;
          for (let i = 0; i < 50; i++) {
            if (root.querySelector('[name="viewer-Filters"] .d4-filter') != null) break;
            await new Promise((r) => setTimeout(r, 100));
          }
          const panelOpen = root.querySelector('[name="viewer-Filters"]') != null;
          const cards = (window as any).__cards(vn) as HTMLElement[];
          const captions = cards.map((e) => (e.querySelector('.d4-filter-column-name')?.textContent ?? '').trim());
          const cardEl = (caption: string) =>
            cards.find((x) => x.querySelector('.d4-filter-column-name')?.textContent?.trim() === caption);
          const raceEl = cardEl('RACE');
          const ageEl = cardEl('AGE');
          const fg = view.getFiltersGroup();
          return {
            panelOpen,
            captions,
            trueCount: view.dataFrame.filter.trueCount,
            raceActive: fg.getStates('RACE', 'categorical')[0]?.active,
            ageActive: fg.getStates('AGE', 'histogram')[0]?.active,
            raceDisabledClass: raceEl ? raceEl.classList.contains('d4-filter-disabled') : null,
            ageDisabledClass: ageEl ? ageEl.classList.contains('d4-filter-disabled') : null,
          };
        }, {id: layoutId, vn: LAYOUT_VIEW, cap: 5000});

        expect(restored.panelOpen).toBe(true);
        expect(restored.captions).toEqual(saved.captions);
        expect(restored.raceActive).toBe(false);
        expect(restored.raceDisabledClass).toBe(true);
        expect(restored.ageActive).toBe(true);
        expect(restored.ageDisabledClass,
          'the AGE card came back painted disabled — the re-apply did not restore per-card state')
          .toBe(false);
        expect(restored.trueCount).toBe(saved.savedTrueCount);
      } finally {
        if (layoutId) {
          await page.evaluate(async (id: string) => {
            try { const s = await grok.dapi.layouts.find(id); await grok.dapi.layouts.delete(s); } catch (_) { /* */ }
          }, layoutId);
        }
      }
    });
  } finally {
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
