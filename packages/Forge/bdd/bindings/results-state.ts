/* What the Train view says about its results, on `.forge-train-results` (ARCHITECTURE.md, "Automation surface"):
   aria-busy while it checks, trains or shows a result, aria-invalid for a problem, and data-forge-run-id, the id of the
   training_run row whose result it shows (absent while it shows a hint); the seed shown and a viewer's reading as text.
   No steps here: the step modules import it. */
import type {Page} from '@playwright/test';
import {type ElementRef, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

export interface ResultsState {busy: string | null; invalid: string | null; run: string | null}

export function resultsState(page: Page): Promise<ResultsState> {
  return page.locator('.forge-train-results').evaluateAll((all) => {
    const e = all.find((x) => x.getClientRects().length > 0);
    return {busy: e?.getAttribute('aria-busy') ?? null, invalid: e?.getAttribute('aria-invalid') ?? null,
      run: e?.getAttribute('data-forge-run-id') ?? null};
  });
}

/** Waits until no Train view on the page has anything pending: no check, no training, no result being shown, whatever
 * its validity. A claim that nothing trained reads after it. */
export async function settled(page: Page, ms = 15000): Promise<void> {
  await expect.poll(() => page.locator('.forge-train-results').evaluateAll((all) =>
    all.every((e) => e.getAttribute('aria-busy') !== 'true')),
  {timeout: pollMs(ms), message: 'the Train view settled (aria-busy not true)'}).toBe(true);
}

/** Where a seed is shown: the Performance panel of the results, or a model's Performance pane in the context panel. */
export const RESULTS_SEED = '.forge-results-metrics .forge-seed';
export const PANE_SEED = '.grok-prop-panel .forge-seed';

/** The seed shown in [where] (the results' Performance panel unless told otherwise). */
export async function shownSeed(page: Page, where = RESULTS_SEED): Promise<string> {
  return (await page.locator(where).filter({visible: true}).first().textContent() ?? '').trim();
}

/** A viewer's reading [name] as text; an absent one as the library's message in brackets, which no value equals. */
export async function readingText(page: Page, target: ElementRef, name: string): Promise<string> {
  const value = await viewers.readingOf(page, target, name);
  return value instanceof viewers.MissingReading ? `(${value.text})` : String(value);
}

/** The run whose result the shown Train view shows, or null. */
export async function shownRun(page: Page): Promise<string | null> {
  return (await resultsState(page)).run;
}
