/* ---
realizes: [filters.cp.save-and-reapply-state]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {cardCount, expectHeaderCounter, expectHeaderCounterQuietNow, trueCount} from '../../helpers/filter-panel';
import {addCardViaPicker} from './column-picker';
import {categoryRowPoint} from './panel-core-ladder-shared';
import {dismissPanelMenu, readCounterTooltipSummary, readSaveOrApplyLeaves, removeProbeState,
  saveStateViaMenu} from './save-and-reapply-shared';

declare const grok: any;

// Steps 1-7 on the local lane. Step 8, whose subject is beer.csv (no local copy), lives in
// save-and-reapply-state-server-spec.ts together with the save it checks against.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const ALWAYS_OFFERED_LEAF = 'Save...';
const probeName = `filter-state-${Date.now()}`;

async function cardCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
      .map((e) => (e.textContent ?? '').trim()));
}

test('Filters — Save and re-apply named filter state', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await openDatagrok(page);
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
      await addCardViaPicker(page, 'RACE');
      expect(await cardCaptions(page),
        'the card that came out of the picker is not captioned RACE — the wrong column got a filter card')
        .toContain('RACE');
      expect(await trueCount(page),
        'adding a card with no category checked already narrowed the table').toBe(total);
      await expectHeaderCounterQuietNow(page,
        'the RACE card filters nothing yet, so the active-filter counter must be hidden or read 0');
    });

    await softStep('Step 2 Real click on the Asian category row → trueCount drops, counter reads 1', async () => {
      const pt = await categoryRowPoint(page, 'RACE', 'Asian');
      expect(pt, 'the RACE card never painted a body tall enough to carry the Asian row').not.toBeNull();
      const clickResult = await page.evaluate(async ({x, y}) => {
        const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
        const race = cards.find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'RACE');
        const overlay = race?.querySelector('[name="overlay"]') as HTMLElement | null;
        if (!overlay) return {before: -1, after: -1};
        const count = () => grok.shell.tv.dataFrame.filter.trueCount;
        const before = count();
        const o = {bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0};
        overlay.dispatchEvent(new MouseEvent('mousedown', o));
        overlay.dispatchEvent(new MouseEvent('mouseup', o));
        overlay.dispatchEvent(new MouseEvent('click', o));
        const after = await (window as any).__moved(count, before, 900);
        return {before, after};
      }, pt!);
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
      await addCardViaPicker(page, 'AGE');
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
      await saveStateViaMenu(page, probeName);
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
        intervals: [30, 60, 120, 250, 500, 1000],
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
  } finally {
    await removeProbeState(page, probeName);
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
