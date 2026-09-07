/* ---
realizes: [filters.cp.expression-and-text-ui]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {trueCount} from '../../helpers/filter-panel';
import {addCardViaPicker} from './column-picker';
import {andOrControl, filterSummaryText, removeAllCards} from './expression-text-shared';

declare const grok: any;

// Scenario 2: the Text filter on beer.csv, which has no local copy.
test.use(specTestOptions);

const beerPath = 'System:DemoFiles/beer.csv';
const fullCount_beer = 118;

async function resetToAromaTextFilter(page: Page): Promise<void> {
  await removeAllCards(page);
  await addCardViaPicker(page, 'Aroma');
  await page.waitForFunction(() => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent ?? '').trim() === 'Aroma');
    if (!card) return false;
    if (!card.querySelector('.d4-text-filter')) return false;
    if (card.querySelector('.d4-update-shadow')) return false;
    return !!card.querySelector('input.d4-search-input') && !!card.querySelector('input[type="range"]');
  }, null, {timeout: 120_000, polling: 60});
}

async function addAromaTerm(page: Page, term: string): Promise<number> {
  return page.evaluate(async (t) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent || '').trim() === 'Aroma') as HTMLElement;
    const input = card.querySelector('.d4-search-input') as HTMLInputElement;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    const before = grok.shell.tv.dataFrame.filter.trueCount;
    input.focus();
    setInp.call(input, t);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    input.dispatchEvent(new KeyboardEvent('keyup', {bubbles: true, key: 'Enter', code: 'Enter', keyCode: 13} as any));
    return (window as any).__moved(() => grok.shell.tv.dataFrame.filter.trueCount, before, 1700);
  }, term);
}

async function setAromaFuzziness(page: Page, value: number): Promise<number> {
  return page.evaluate(async (val) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent || '').trim() === 'Aroma') as HTMLElement;
    const range = card.querySelector('input[type="range"]') as HTMLInputElement;
    const setInp = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    // The text filter matches asynchronously, so a settle that gives up as soon as the frame has
    // raised no filter pass reads the pre-slider count: the whole budget is the wait here.
    const before = grok.shell.tv.dataFrame.filter.trueCount;
    range.focus();
    setInp.call(range, String(val));
    return (window as any).__filtered(() => {
      range.dispatchEvent(new Event('input', {bubbles: true}));
      range.dispatchEvent(new Event('change', {bubbles: true}));
    }, 1800, before);
  }, value);
}

test('Filter Panel — Text filter driven through its own UI', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.openTable(page, {path: beerPath, withFilterPanel: true, semTypeTimeoutMs: 1000});

  let trueCount_or_beer = 0;
  let trueCount_and_beer = 0;

  await softStep('Scenario 2 Step 13 — Aroma text filter: a term drops the beer row count', async () => {
    await resetToAromaTextFilter(page);
    expect(await trueCount(page), 'baseline is the full beer row set').toBe(fullCount_beer);
    await setAromaFuzziness(page, 0);
    const afterMalt = await addAromaTerm(page, 'malt');
    expect(afterMalt, 'typing a term and pressing Enter drops the beer row count').toBeLessThan(fullCount_beer);
    expect(afterMalt).toBeGreaterThan(0);
    expect(await filterSummaryText(page), 'the summary reflects the active Aroma term').toContain('malt');

    const content = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('Aroma');
      let passingWithTerm = 0, passingWithoutTerm = 0, excludedWithTerm = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const has = String(col.get(i) ?? '').toLowerCase().includes('malt');
        if (df.filter.get(i)) has ? passingWithTerm++ : passingWithoutTerm++;
        else if (has) excludedWithTerm++;
      }
      return {passingWithTerm, passingWithoutTerm, excludedWithTerm};
    });
    expect(content.passingWithTerm, 'no row containing the term survived — the term is not the criterion')
      .toBeGreaterThan(0);
    expect(content.passingWithoutTerm, 'rows whose Aroma does not contain the term passed the filter')
      .toBe(0);
    expect(content.excludedWithTerm, 'rows whose Aroma contains the term were filtered out').toBe(0);
  });

  await softStep('Scenario 2 Step 14 — two terms: AND is strictly stricter than OR', async () => {
    await addAromaTerm(page, 'hop');
    expect(await andOrControl(page, 'Aroma', false), 'the two-term card starts in OR mode').toBe('OR');
    trueCount_or_beer = await trueCount(page);
    const afterToggle = await andOrControl(page, 'Aroma', true);
    expect(afterToggle, 'the header toggle switches OR to AND').toBe('AND');
    trueCount_and_beer = await trueCount(page);
    expect(trueCount_and_beer, 'AND is strictly stricter than OR for the same two terms — the control is wired')
      .toBeLessThan(trueCount_or_beer);
  });

  await softStep('Scenario 2 Step 15 — fuzziness 0 yields zero matches; raising the slider raises the count above 0 and grows it', async () => {
    await resetToAromaTextFilter(page);
    expect(await trueCount(page),
      'the fresh Aroma card must be carrying no criterion — without a proven full beer row set here, a leftover filter already sitting at zero rows would satisfy the near-miss check below on its own')
      .toBe(fullCount_beer);
    const atReset = await setAromaFuzziness(page, 0);
    expect(atReset, 'pinning fuzziness to 0 filters nothing on its own').toBe(fullCount_beer);
    const atZero = await addAromaTerm(page, 'maltx');
    expect(atZero, 'at fuzziness 0 a non-matching near-miss term yields zero matches').toBe(0);
    const atMid = await setAromaFuzziness(page, 0.5);
    const atHigh = await setAromaFuzziness(page, 0.8);
    expect(atMid, 'raising fuzziness recovers matches — the count rises above 0').toBeGreaterThan(0);
    expect(atHigh, 'raising fuzziness further grows the matched set').toBeGreaterThan(atMid);
    expect(atHigh).toBeLessThanOrEqual(fullCount_beer);
  });

  await softStep('Teardown', async () => {
    await v.closeAllAndWait(page);
  });

  v.finishSpec();
});
