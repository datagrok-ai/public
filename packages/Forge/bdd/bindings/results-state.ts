/* What the Train view says about its results, on `.forge-train-results` (ARCHITECTURE.md, "Automation surface"):
   aria-busy while it checks, trains or shows a result, aria-invalid for a problem, and data-forge-run-id, the id of the
   training_run row whose result it shows (absent while it shows a hint). No steps here: the step modules import it. */
import type {Page} from '@playwright/test';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

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

/** The seed in the shown Performance panel. */
export async function shownSeed(page: Page): Promise<string> {
  return (await page.locator('.forge-results-metrics .forge-seed').filter({visible: true}).first().textContent() ?? '').trim();
}

/** The run whose result the shown Train view shows, or null. */
export async function shownRun(page: Page): Promise<string | null> {
  return (await resultsState(page)).run;
}
