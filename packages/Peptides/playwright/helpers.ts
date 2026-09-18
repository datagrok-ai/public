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
 * The wait ends early once `startAnalysis` closes its "Loading SAR..." progress indicator:
 * whatever is attached then is all that will be. Waiting on a fixed budget instead read the
 * viewers mid-wave on a slow runner — Logo Summary Table takes ~40 s on dev and longer on a
 * two-core CI agent, which is what failed the SAR specs on GitHub Actions.
 *
 * A timeout is NOT thrown: the caller's `expect` names the viewer that is missing, which
 * is a far better failure message than a bare waitForFunction timeout.
 */
export async function waitForViewers(page: Page, expected: string[], timeoutMs = 300_000): Promise<string[]> {
  const started = Date.now();
  while (Date.now() - started < timeoutMs) {
    // TEMPORARY instrumentation (revert before merge): the SAR specs fail on GitHub Actions
    // with Logo Summary Table never attached, while dev attaches it at 39-44 s. This says
    // where the analysis stops on a two-core agent.
    const probe = await page.evaluate((want: string[]) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv?.dataFrame?.temp['peptidesModel'];
      const types = Array.from(tv?.viewers ?? []).map((v) => v.type);
      const bars = Array.from(document.querySelectorAll('.d4-task-bar, .d4-progress'))
        .map((e) => (e.textContent ?? '').trim()).filter((t) => t.length > 0);
      return {
        done: want.every((t) => types.includes(t)),
        types, bars,
        mclCols: (model?._mclCols ?? []) as string[],
        columns: (tv?.dataFrame?.columns?.names() ?? []).filter((n: string) => /cluster|embed/i.test(n)),
      };
    }, expected);
    const elapsed = Math.round((Date.now() - started) / 1000);
    console.log(`[sar-probe] ${elapsed}s viewers=${JSON.stringify(probe.types)} bars=${JSON.stringify(probe.bars)} ` +
      `mclCols=${JSON.stringify(probe.mclCols)} clusterCols=${JSON.stringify(probe.columns)}`);
    if (probe.done)
      break;
    if (elapsed > 5 && !probe.bars.some((t) => t.includes('Loading SAR')))
      break;
    await page.waitForTimeout(15_000);
  }
  return attachedViewers(page);
}
