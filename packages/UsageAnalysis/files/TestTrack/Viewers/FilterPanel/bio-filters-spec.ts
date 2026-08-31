/* ---
realizes: [filters.cp.chem-and-bio-filters]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const BIO_PATH = 'System:DemoFiles/bio/peptides.csv';
const BIO_SUBSEQ = 'T-T-Y-K-N-Y-V';
const PICKER = '.d4-dialog[name="dialog-Select-columns..."]';
const BIO_REGISTERED_WINDOW_MS = 20_000;
const BIO_READY_WINDOW_MS = 60_000;
const SEQ_HELPER_NOT_INITIALIZED = 'seqHelper is not initialized';

interface BioFactoryProbe {
  found: boolean;
  applied: boolean;
  error: string;
  fatal: boolean;
  foundMs: number;
  appliedMs: number;
  waitedMs: number;
  notInitRetries: number;
  falsyRetries: number;
}


async function probeBioFactory(page: Page): Promise<BioFactoryProbe> {
  return page.evaluate(async ({windowMs, registeredWindowMs, retryablePattern}) => {
    const retryable = new RegExp(retryablePattern, 'i');
    const t0 = Date.now();
    const since = () => Date.now() - t0;
    let found = false;
    let foundMs = -1;
    let error = '';
    let notInitRetries = 0;
    let falsyRetries = 0;
    while (since() < windowMs) {
      const fns = DG.Func.find({name: 'bioSubstructureFilter'}) || [];
      const f = fns.find((x: any) => x.package?.name === 'Bio');
      if (f) {
        if (!found) {
          found = true;
          foundMs = since();
        }
        try {
          if (await f.apply())
            return {found: true, applied: true, error: '', fatal: false,
              foundMs, appliedMs: since(), waitedMs: since(), notInitRetries, falsyRetries};
          error = 'apply() resolved to a falsy widget';
          falsyRetries++;
        } catch (e: any) {
          error = String(e?.message ?? e);
          if (!retryable.test(error))
            return {found: true, applied: false, error, fatal: true,
              foundMs, appliedMs: -1, waitedMs: since(), notInitRetries, falsyRetries};
          notInitRetries++;
        }
      }
      await new Promise((r) => setTimeout(r, 500));
      if (!found && since() >= registeredWindowMs) break;
    }
    return {found, applied: false, error, fatal: false,
      foundMs, appliedMs: -1, waitedMs: since(), notInitRetries, falsyRetries};
  }, {windowMs: BIO_READY_WINDOW_MS, registeredWindowMs: BIO_REGISTERED_WINDOW_MS,
    retryablePattern: SEQ_HELPER_NOT_INITIALIZED});
}

async function macromoleculeColumn(page: Page): Promise<string | null> {
  return page.evaluate(async () => {
    for (let i = 0; i < 40; i++) {
      const df = grok.shell.tv.dataFrame;
      for (let j = 0; j < df.columns.length; j++) {
        const c = df.columns.byIndex(j);
        if (c.semType === 'Macromolecule') return c.name as string;
      }
      await new Promise((r) => setTimeout(r, 500));
    }
    return null;
  });
}

async function bioCardCount(page: Page, column: string): Promise<number> {
  return page.evaluate((col) => {
    let cards = 0;
    for (const card of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
      const caption = (card.querySelector('.d4-filter-column-name')?.textContent ?? '').trim();
      if (caption === col && card.querySelector('input[placeholder="Substructure"]')) cards++;
    }
    return cards;
  }, column);
}

async function waitForPanelSettled(page: Page): Promise<number> {
  return page.evaluate(async () => {
    const count = () => document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length;
    let stable = 0;
    let last = -1;
    for (let i = 0; i < 60; i++) {
      const now = count();
      stable = now === last ? stable + 1 : 0;
      last = now;
      if (stable >= 4) return now;
      await new Promise((r) => setTimeout(r, 500));
    }
    return last;
  });
}

async function filterGroupSize(page: Page): Promise<number> {
  return page.evaluate(() => [...grok.shell.tv.getFiltersGroup().filters].length);
}

async function removeSubstructureCards(page: Page): Promise<number> {
  return page.evaluate(() => {
    let clicks = 0;
    for (const card of document.querySelectorAll('[name="viewer-Filters"] .d4-filter')) {
      if (!card.querySelector('input[placeholder="Substructure"]')) continue;
      const x = card.querySelector('[name="icon-times"]') as HTMLElement | null;
      if (!x) throw new Error('a substructure filter card carries no [name="icon-times"] remove icon');
      x.click();
      clicks++;
    }
    return clicks;
  });
}

async function confirmColumnPicker(page: Page): Promise<{before: number; after: number}> {
  const dialog = page.locator(PICKER);
  await dialog.waitFor({timeout: 20_000});
  const readChecked = async () => dialog.evaluate((d) => {
    const label = [...d.querySelectorAll('label')]
      .find((l) => /^\s*\d+\s+checked/.test(l.textContent ?? ''));
    const m = (label?.textContent ?? '').trim().match(/^(\d+)\s+checked/);
    return m ? parseInt(m[1], 10) : -1;
  });
  const before = await readChecked();
  await dialog.locator('[name="label-All"]').click();
  await expect.poll(readChecked, {
    timeout: 15_000,
    message: 'the "All" link of the Select columns... picker never moved its "<n> checked" counter off ' +
      'zero (-1 means the counter label itself is gone), so nothing was offered to check and OK would ' +
      'create no card',
  }).toBeGreaterThan(0);
  const after = await readChecked();
  await dialog.locator('[name="button-OK"]').click();
  await dialog.waitFor({state: 'detached', timeout: 20_000});
  return {before, after};
}

async function trueCount(page: Page): Promise<number> {
  const c = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
  expect(typeof c, 'trueCount must be a number, not null/undefined').toBe('number');
  return c;
}

test('Filter panel — Bio package filters', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  try {
    await v.openTable(page, {path: BIO_PATH, semType: 'Macromolecule'});
    const bioFullCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(bioFullCount, 'peptides.csv must open with rows').toBeGreaterThan(0);

    const factory = await probeBioFactory(page);
    console.log(`Bio full count: ${bioFullCount}; Bio factory registered: ${factory.found} ` +
      `at ${factory.foundMs} ms; applicable: ${factory.applied} at ${factory.appliedMs} ms; ` +
      `waited ${factory.waitedMs} ms of ${BIO_READY_WINDOW_MS} ms; ` +
      `retryable not-initialized sightings: ${factory.notInitRetries}; ` +
      `falsy-widget retries: ${factory.falsyRetries}` +
      `${factory.error ? `; last error: ${factory.error}` : ''}`);
    test.skip(!factory.found,
      `no Bio-owned bioSubstructureFilter function is registered within ${BIO_REGISTERED_WINDOW_MS / 1000}s ` +
      '— the Bio package is not installed on this build, so the Bio filter scenario has no subject');
    expect(factory.applied,
      `Bio:bioSubstructureFilter IS registered (seen at ${factory.foundMs} ms) but never became applicable, ` +
      `so the card the panel menu is about to build cannot be created. Waited ${factory.waitedMs} ms of the ` +
      `${BIO_READY_WINDOW_MS} ms readiness window, over which "${SEQ_HELPER_NOT_INITIALIZED}" was tolerated ` +
      `as a startup race ${factory.notInitRetries} time(s) (and a falsy widget ${factory.falsyRetries} ` +
      'time(s)); the last error was ' +
      `${factory.error || 'no error was raised'}. A high retry count that ran the window out means the ` +
      'window is too short for this build; a fatal error on the first sighting means the package itself ' +
      `is broken (fatal: ${factory.fatal}).`)
      .toBe(true);

    const macroColumn = await macromoleculeColumn(page);
    expect(macroColumn,
      'no column of peptides.csv reports semType Macromolecule, so the Bio substructure card has no ' +
      'column to bind to and the picker below would offer nothing')
      .not.toBeNull();
    console.log(`Macromolecule column: ${macroColumn}`);

    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"]').first().waitFor({timeout: 30_000});
    await page.waitForFunction(() => {
      const fv = [...(grok.shell.tv?.viewers ?? [])].find((x: any) => x.type === 'Filters');
      return !!(fv && (fv as any).props);
    }, null, {timeout: 30_000, polling: 200});
    console.log(`Filter Panel settled at ${await waitForPanelSettled(page)} card(s)`);

    await softStep('Scenario 1, Step 3: Enter subsequence in Bio substructure filter; trueCount drops below full', async () => {
      const removed = await removeSubstructureCards(page);
      await expect.poll(() => bioCardCount(page, macroColumn!), {
        timeout: 15_000,
        intervals: [400, 800, 1500],
        message: `a substructure card bound to ${macroColumn} survived ${removed} remove-icon click(s), ` +
          'so the Add Filter > Bio Substructure Filter... drive below could be satisfied by a card it ' +
          'did not create',
      }).toBe(0);
      const groupSizeBefore = await filterGroupSize(page);

      await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Bio Substructure Filter...');
      const picker = await confirmColumnPicker(page);
      expect(picker.before,
        `the Select columns... picker opened with ${picker.before} column(s) already checked (-1 means ` +
        'its counter label is gone), so the "All" click below proves nothing about what it offered')
        .toBe(0);
      expect(picker.after,
        'the Select columns... picker reported no column checked after "All", so OK confirmed an empty ' +
        'selection and the Macromolecule column was never offered')
        .toBeGreaterThan(0);

      await expect.poll(() => bioCardCount(page, macroColumn!), {
        timeout: 45_000,
        intervals: [500, 1000, 2000],
        message: 'Add Filter > Bio Substructure Filter... did not take the panel from zero to exactly one ' +
          `Bio substructure card bound to ${macroColumn} — either the menu leaf created nothing, or the ` +
          'card it created is bound to a different column, or the leaf fired twice',
      }).toBe(1);
      expect(await filterGroupSize(page),
        `the Filters group held ${groupSizeBefore} filter(s) before the Bio leaf was driven and does not ` +
        'hold exactly one more after it — the card on screen is not backed by a new filter in the group')
        .toBe(groupSizeBefore + 1);

      const cardSel = `[name="viewer-Filters"] .d4-filter:has(.d4-filter-column-name:text-is("${macroColumn}"))`;
      const subLocator = page.locator(`${cardSel} input[placeholder="Substructure"]`);
      await subLocator.waitFor({timeout: 45_000});
      await subLocator.fill(BIO_SUBSEQ);
      await subLocator.press('Enter');
      await expect.poll(() => subLocator.inputValue(), {
        timeout: 15_000,
        message: 'the Substructure input does not read back the subsequence that was typed into it — ' +
          'the card re-rendered and the text landed on an element that is no longer on screen',
      }).toBe(BIO_SUBSEQ);
      await page.waitForFunction(
        (full) => grok.shell.tv.dataFrame.filter.trueCount < full, bioFullCount, {timeout: 30_000});
      const active = await trueCount(page);
      expect(active).toBeLessThan(bioFullCount);
      expect(active, 'the probe subsequence must keep at least one peptide row').toBeGreaterThan(0);
      console.log(`Bio subsequence active trueCount: ${active}`);
    });

    await softStep('Scenario 1, Step 4: Clear the Bio substructure filter; all peptide rows return', async () => {
      const cardSel = `[name="viewer-Filters"] .d4-filter:has(.d4-filter-column-name:text-is("${macroColumn}"))`;
      const subLocator = page.locator(`${cardSel} input[placeholder="Substructure"]`);
      expect(await trueCount(page), 'the filter must still be narrowing before it is cleared').toBeLessThan(bioFullCount);
      await expect.poll(() => subLocator.inputValue(), {
        timeout: 15_000,
        intervals: [300, 600, 1000],
        message: 'the Substructure input no longer holds the subsequence Step 3 typed, so clearing it ' +
          'would be clearing something the product had already dropped on its own',
      }).toBe(BIO_SUBSEQ);
      await subLocator.fill('');
      await subLocator.press('Enter');
      await expect.poll(() => subLocator.inputValue(), {
        timeout: 15_000,
        message: 'the Substructure input did not stay empty after it was cleared — the card re-rendered ' +
          'and restored the previous text, so the clear never reached the live element',
      }).toBe('');
      await page.waitForFunction(
        (full) => grok.shell.tv.dataFrame.filter.trueCount === full, bioFullCount, {timeout: 30_000});
      expect(await trueCount(page)).toBe(bioFullCount);
    });
  } finally {
    await v.cleanupShell(page);
  }

  v.finishSpec('Bio scenario step failures');
});
