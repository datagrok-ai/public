/**
 * Playwright helpers shared by Bio/* spec files.
 *
 * Each helper is a verbatim extraction of a block previously pasted into the
 * Bio specs — same selectors, same mechanism, same sleeps. Imported as:
 * `import * as bio from '../helpers/bio';`.
 */

import {Page} from '@playwright/test';

/**
 * Launch a Bio > Analyze leaf via the top menu. Verbatim equivalent of the
 * block in the Bio specs: click `[name="div-Bio"]` → 400ms → mouseover
 * `[name="div-Bio---Analyze"]` (hover-not-click surfaces the leaves) → 300ms →
 * click the leaf. `leaf` is the leaf element's `name` (e.g.
 * `div-Bio---Analyze---Sequence-Space...`).
 */
export async function openBioAnalyze(page: Page, leaf: string): Promise<void> {
  await page.evaluate(async (leafSel) => {
    (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
    await new Promise((r) => setTimeout(r, 400));
    document.querySelector('[name="div-Bio---Analyze"]')!
      .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 300));
    (document.querySelector('[name="' + leafSel + '"]') as HTMLElement).click();
  }, leaf);
}
