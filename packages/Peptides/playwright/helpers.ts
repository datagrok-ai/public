import {Page} from '@playwright/test';

/** Viewer types attached to the TableView that owns the PeptidesModel. */
export function attachedViewers(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
    return Array.from(tv.viewers).map((v) => v.type);
  });
}

/**
 * Wait until every viewer in `expected` is attached, then return the whole set.
 *
 * SAR attaches its viewers in waves: the model appears first, MCL clustering finishes
 * later, and Logo Summary Table only exists once clustering produced a clusters column.
 * A fixed settle-and-read reports whatever happened to be mounted at that instant, which
 * is how the Logo Summary Table steps used to pass while asserting nothing.
 *
 * A timeout is NOT thrown: the caller's `expect` names the viewer that is missing, which
 * is a far better failure message than a bare waitForFunction timeout.
 */
export async function waitForViewers(page: Page, expected: string[], timeoutMs = 120_000): Promise<string[]> {
  try {
    await page.waitForFunction((want: string[]) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const types = Array.from(tv.viewers).map((v) => v.type);
      return want.every((t) => types.includes(t));
    }, expected, {timeout: timeoutMs});
  } catch (_) {
  }
  return attachedViewers(page);
}
