/* ---
realizes: [filters.cp.save-and-reapply-state]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addCardViaColumnSelector, cardCount, expectHeaderCounter, expectHeaderCounterQuietNow,
  HEADER_COUNTER, trueCount} from '../../helpers/filter-panel';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const beerPath = 'System:DemoFiles/beer.csv';
const AMBIENT_CONSOLE_ERROR = 'Permissions policy violation: compute-pressure';
const ALWAYS_OFFERED_LEAF = 'Save...';
const probeName = `filter-state-${Date.now()}`;

async function readCounterTooltipSummary(page: Page): Promise<Record<string, string> | null> {
  return page.evaluate(async (sel) => {
    const counterEl = document.querySelector(sel);
    if (!counterEl) return null;
    const r = counterEl.getBoundingClientRect();
    const at = {clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, bubbles: true};
    counterEl.dispatchEvent(new MouseEvent('mouseenter', at));
    counterEl.dispatchEvent(new MouseEvent('mousemove', at));
    counterEl.dispatchEvent(new MouseEvent('mouseover', at));
    for (let i = 0; i < 20; i++) {
      if (document.querySelector('.d4-tooltip table')) break;
      await new Promise((res) => setTimeout(res, 100));
    }
    const table = document.querySelector('.d4-tooltip table');
    if (!table) return null;
    const summary: Record<string, string> = {};
    for (const row of Array.from(table.querySelectorAll('tr'))) {
      const cells = Array.from(row.querySelectorAll('td,th'));
      if (cells.length < 2) continue;
      const caption = cells[0].textContent?.trim() ?? '';
      if (caption.length > 0) summary[caption] = cells[1].textContent?.trim() ?? '';
    }
    return Object.keys(summary).length > 0 ? summary : null;
  }, HEADER_COUNTER);
}

async function cardCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
      .map((e) => (e.textContent ?? '').trim()));
}

async function readSaveOrApplyLeaves(page: Page): Promise<string[]> {
  await page.evaluate(() => {
    const root = grok.shell.tv.root as HTMLElement;
    const filtersViewer = root.querySelector('[name="viewer-Filters"]');
    const titlebar = filtersViewer?.closest('.panel-base')?.querySelector('.panel-titlebar');
    const hamburger = titlebar?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement | null;
    if (!hamburger) {
      throw new Error('no [name="icon-font-icon-menu"] in the Filter Panel title bar of view ' +
        `"${grok.shell.tv.name}"`);
    }
    hamburger.click();
  });
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  return page.evaluate(async () => {
    const popup = () => Array.from(document.querySelectorAll('.d4-menu-popup')).pop() ?? null;
    const labelsOf = (root: Element | null) =>
      Array.from(root?.querySelectorAll('.d4-menu-item-label') ?? [])
        .map((l) => (l.textContent ?? '').trim());
    const findGroup = () => Array.from(popup()?.children ?? [])
      .filter((c) => c.classList.contains('d4-menu-item'))
      .find((it) => it.querySelector(':scope > .d4-menu-item-label')?.textContent?.trim() ===
        'Save or Apply') ?? null;
    const groupDeadline = Date.now() + 5000;
    let group = findGroup();
    while (!group && Date.now() < groupDeadline) {
      await new Promise((r) => setTimeout(r, 100));
      group = findGroup();
    }
    if (!group) {
      throw new Error('panel menu: group "Save or Apply" not found in the panel popup; visible: ' +
        labelsOf(popup()).join(' | '));
    }
    const rect = group.getBoundingClientRect();
    const at = {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true};
    for (const type of ['mouseover', 'mouseenter', 'mousemove'])
      group.dispatchEvent(new MouseEvent(type, at));
    const leaves = () => Array.from(group!.querySelectorAll('.d4-menu-item'))
      .map((c) => c.querySelector(':scope > .d4-menu-item-label')?.textContent?.trim() ?? '')
      .filter((s) => s.length > 0);
    const leafDeadline = Date.now() + 5000;
    let out = leaves();
    while (out.length === 0 && Date.now() < leafDeadline) {
      await new Promise((r) => setTimeout(r, 100));
      out = leaves();
    }
    return out;
  });
}

async function visibleMenuPopups(page: Page): Promise<number> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-menu-popup')).filter((p) => {
      const cs = window.getComputedStyle(p as HTMLElement);
      const r = (p as HTMLElement).getBoundingClientRect();
      return cs.display !== 'none' && cs.visibility !== 'hidden' && r.width > 0 && r.height > 0;
    }).length);
}

async function dismissPanelMenu(page: Page): Promise<void> {
  for (let attempt = 0; attempt < 3; attempt++) {
    await page.keyboard.press('Escape');
    const deadline = Date.now() + 3000;
    while (Date.now() < deadline) {
      if (await visibleMenuPopups(page) === 0) return;
      await page.waitForTimeout(200);
    }
  }
  expect(await visibleMenuPopups(page),
    'the panel menu popup stayed on screen and would intercept the next gesture').toBe(0);
}

test('Filters — Save and re-apply named filter state', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  let total = -1;
  let trueCountSaved = -1;
  let summarySaved: Record<string, string> | null = null;

  try {
    await softStep('Setup demog opens with a Filter Panel holding no cards', async () => {
      total = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(typeof total, 'demog did not report a row count — every later comparison would be vacuous')
        .toBe('number');
      expect(total, 'demog opened empty, so "narrower than the full table" could not discriminate')
        .toBeGreaterThan(0);
      await v.resetFilters(page);
      const headers = await page.locator('[name="viewer-Filters"] .d4-filter-group-header').count();
      expect(headers, 'the Filter Panel itself is not on screen after the reset, so "no cards left" ' +
        'would be satisfied by the panel being gone rather than by an empty panel').toBe(1);
      expect(await cardCount(page), 'the panel still holds filter cards after the reset, so the ' +
        'counter cannot be attributed to the two cards this scenario configures').toBe(0);
      expect(await trueCount(page), 'the reset left rows filtered out of demog').toBe(total);
    });

    await softStep('Step 1 Add RACE categorical filter card', async () => {
      await addCardViaColumnSelector(page, 'RACE');
      expect(await cardCaptions(page),
        'the card that came out of the picker is not captioned RACE — the wrong column got a filter card')
        .toContain('RACE');
      expect(await trueCount(page),
        'adding a card with no category checked already narrowed the table').toBe(total);
      await expectHeaderCounterQuietNow(page,
        'the RACE card filters nothing yet, so the active-filter counter must be hidden or read 0');
    });

    await softStep('Step 2 Real click on the Asian category row → trueCount drops, counter reads 1', async () => {
      const clickResult = await page.evaluate(async () => {
        const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
        const race = cards.find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE');
        if (!race) return {before: -1, after: -1};
        const canvas = race.querySelector('[name="canvas"]') as HTMLCanvasElement | null;
        const overlay = race.querySelector('[name="overlay"]') as HTMLElement | null;
        if (!canvas || !overlay) return {before: -1, after: -1};
        const cats: string[] = grok.shell.tv.dataFrame.col('RACE').categories;
        const rowIndex = cats.indexOf('Asian');
        const rect = canvas.getBoundingClientRect();
        const rowHeight = rect.height / cats.length;
        const clientX = rect.left + 10;
        const clientY = rect.top + rowIndex * rowHeight + rowHeight / 2;
        const before = grok.shell.tv.dataFrame.filter.trueCount;
        const o = {bubbles: true, cancelable: true, view: window, clientX, clientY, button: 0};
        overlay.dispatchEvent(new MouseEvent('mousedown', o));
        overlay.dispatchEvent(new MouseEvent('mouseup', o));
        overlay.dispatchEvent(new MouseEvent('click', o));
        await new Promise((r) => setTimeout(r, 900));
        const after = grok.shell.tv.dataFrame.filter.trueCount;
        return {before, after};
      });
      expect(clickResult.before,
        'the RACE card, its canvas or its overlay was not on screen, so no click was delivered')
        .toBeGreaterThanOrEqual(0);
      expect(clickResult.after,
        'the click on the Asian row left the filtered row count unchanged — the gesture missed the row')
        .not.toBe(clickResult.before);
      const count = await trueCount(page);
      expect(count, 'the Asian selection did not narrow the table below its full row count')
        .toBeLessThan(total);
      expect(count, 'the Asian selection filtered every row out').toBeGreaterThan(0);
      await expectHeaderCounter(page, '1',
        'one card started filtering, so the header active-filter counter must read 1');
    });

    await softStep('Step 3 Add AGE histogram 30-60 → counter reads 2, count below full', async () => {
      await addCardViaColumnSelector(page, 'AGE');
      const count = await v.applyNumericFilter(page, 'AGE', 30, 60);
      await expectHeaderCounter(page, '2',
        'RACE and AGE are both filtering, so the header active-filter counter must read 2');
      expect(count, 'the two configured filters did not narrow demog below its full row count')
        .toBeLessThan(total);
      expect(count, 'the two configured filters left no rows, so the round-trip value carries no signal')
        .toBeGreaterThan(0);
      trueCountSaved = count;

      summarySaved = await readCounterTooltipSummary(page);
      expect(summarySaved, 'the counter tooltip rendered no summary table, so no criteria baseline exists')
        .not.toBeNull();
      expect(Object.keys(summarySaved!).sort(),
        'the counter tooltip does not summarise exactly the two configured columns')
        .toEqual(['AGE', 'RACE']);
      expect(summarySaved!['RACE'], 'the RACE criteria summary does not name the Asian category')
        .toContain('Asian');
      expect(summarySaved!['AGE'],
        'the AGE criteria summary is not a rendered numeric range — there is nothing to compare at Step 7')
        .toMatch(/^\[.+,.+\]$/);
    });

    await softStep('Step 4 Save state via hamburger → Save or Apply → Save...', async () => {
      const consoleErrors: string[] = [];
      const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
        if (msg.type() === 'error') consoleErrors.push(msg.text());
      };
      page.on('console', onConsole);
      const errorSet = () => new Set(consoleErrors.filter((t) => !t.includes(AMBIENT_CONSOLE_ERROR)));
      const errorsBefore = errorSet();

      await v.drivePanelMenuLeaf(page, 'Filters', 'Save or Apply', 'Save...');

      const dlg = page.locator('.d4-dialog[name="dialog-Save-filter-state"]');
      await dlg.waitFor({timeout: 10_000});
      const nameInput = dlg.locator('input[name="input-Name"]');
      await nameInput.click();
      await nameInput.fill('');
      await nameInput.pressSequentially(probeName, {delay: 15});
      await dlg.locator('[name="button-OK"]').click();
      await expect.poll(async () => dlg.count(), {
        timeout: 10_000,
        message: 'the Save filter state dialog did not close after OK',
      }).toBe(0);

      const stored = await page.evaluate((name) => {
        const raw = window.localStorage.getItem('filter-states');
        if (raw === null) return {parsed: false, has: false};
        const states = JSON.parse(raw);
        return {parsed: true, has: Object.prototype.hasOwnProperty.call(states, name)};
      }, probeName);
      expect(stored.parsed,
        'the "filter-states" localStorage entry is absent after the save — nothing was stored').toBe(true);
      expect(stored.has,
        'the probe name is absent from the "filter-states" localStorage entry — the state was not saved')
        .toBe(true);

      page.off('console', onConsole);
      const newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
      expect(newErrors, `the save produced new console error texts: ${newErrors.join(' | ')}`).toEqual([]);
    });

    await softStep('Step 5 Perturb filters → trueCount differs from saved', async () => {
      await v.applyNumericFilter(page, 'AGE', 0, 200);
      const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', ['Black']);
      expect(filteredCount,
        'the perturbation left the filtered row count where the saved state had it, so a re-apply that ' +
        'did nothing would still look like a successful round-trip')
        .not.toBe(trueCountSaved);
    });

    await softStep('Step 6 Saved probe name offered under Save or Apply', async () => {
      const leaves = await readSaveOrApplyLeaves(page);
      expect(leaves,
        'the Save or Apply submenu enumerated nothing on demog, so neither presence nor absence of a ' +
        'named state can be read off it')
        .toContain(ALWAYS_OFFERED_LEAF);
      expect(leaves,
        'the saved probe name is not offered under Save or Apply on demog, whose columns it was saved from')
        .toContain(probeName);
      await dismissPanelMenu(page);
    });

    await softStep('Step 7 Re-apply via menu → round-trip restored, no re-entrancy error', async () => {
      const consoleErrors: string[] = [];
      const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
        if (msg.type() === 'error') consoleErrors.push(msg.text());
      };
      page.on('console', onConsole);

      await v.drivePanelMenuLeaf(page, 'Filters', 'Save or Apply', probeName);

      await expect.poll(async () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount), {
        timeout: 15_000,
        message: 'the re-applied state did not restore the filtered row count recorded before the perturbation',
      }).toBe(trueCountSaved);
      await expectHeaderCounter(page, '2',
        'the re-applied state puts two cards back into filtering, so the counter must read 2');

      const summary = await readCounterTooltipSummary(page);
      expect(summary, 'the counter tooltip rendered no summary table after the re-apply').not.toBeNull();
      expect(Object.keys(summary!).sort(),
        'the re-applied state summarises a different set of columns than the saved one')
        .toEqual(['AGE', 'RACE']);
      expect(summary, 'the re-applied criteria differ from the ones saved at Step 3').toEqual(summarySaved);

      page.off('console', onConsole);
      const reentrancy = consoleErrors.filter((t) => /Cannot fire new event/i.test(t));
      expect(reentrancy,
        'a "Cannot fire new event" re-entrancy error surfaced while re-applying through the real menu ' +
        '(GROK-20386)')
        .toEqual([]);
    });

    await softStep('Step 8 Probe name absent under Save or Apply on beer table', async () => {
      await page.evaluate(async (path) => {
        const beer = await grok.dapi.files.readCsv(path);
        const tv2 = grok.shell.addTableView(beer);
        await new Promise((resolve) => {
          const sub = beer.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
          setTimeout(resolve, 3000);
        });
        tv2.getFiltersGroup();
      }, beerPath);
      await expect.poll(async () => page.evaluate(() =>
        (grok.shell.tv.root as HTMLElement)
          .querySelectorAll('[name="viewer-Filters"] .d4-filter-group-header').length), {
        timeout: 20_000,
        intervals: [500, 1000, 2000],
        message: 'the beer table view never grew a Filter Panel to read the Save or Apply submenu from',
      }).toBe(1);
      const leaves = await readSaveOrApplyLeaves(page);
      expect(leaves,
        'the Save or Apply submenu enumerated nothing on beer.csv, so the absence of the probe name ' +
        'below would be satisfied by a submenu that simply failed to populate')
        .toContain(ALWAYS_OFFERED_LEAF);
      expect(leaves,
        'the demog-typed state is offered on beer.csv, whose columns do not match it')
        .not.toContain(probeName);
      await dismissPanelMenu(page);
    });
  } finally {
    await page.evaluate((name) => {
      try {
        const states = JSON.parse(window.localStorage.getItem('filter-states') ?? '{}');
        if (Object.prototype.hasOwnProperty.call(states, name)) {
          delete states[name];
          window.localStorage.setItem('filter-states', JSON.stringify(states));
        }
      } catch (_) {}
    }, probeName);
    await softStep('Teardown probe key removed from localStorage', async () => {
      const leaked = await page.evaluate((name) => {
        const raw = window.localStorage.getItem('filter-states');
        if (raw === null) return false;
        const states = JSON.parse(raw);
        return Object.prototype.hasOwnProperty.call(states, name);
      }, probeName);
      expect(leaked, 'the probe named state is still in localStorage and would bleed into later runs')
        .toBe(false);
    });
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
