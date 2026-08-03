/**
 * Playwright helpers shared by Chem/* spec files.
 *
 * Each helper is a verbatim extraction of a block previously pasted into the
 * Chem specs — same selectors, same mechanism, same sleeps. Imported as:
 * `import * as chem from '../helpers/chem';` (or `'../../helpers/chem'` from
 * Chem/Advanced/*).
 */

import {Page} from '@playwright/test';

/**
 * Open a Chem top-menu item by visible label. Verbatim equivalent of the
 * menu-navigation block in the Chem specs: dispatch a click on
 * `[name="div-Chem"]` → `delayMs` (submenu render) → find the
 * `.d4-menu-item-label` whose trimmed text === `label` → dispatch a click on
 * its `.closest('.d4-menu-item')`. The caller keeps its own post-open wait
 * (`.d4-dialog` / viewer probe) and assertions. `delayMs` defaults to 600 (the
 * most common site); pass each site's exact delay (600/800/…) via opts.
 */
export async function openChemMenuItem(
  page: Page, label: string, opts?: {delayMs?: number},
): Promise<void> {
  const delayMs = opts?.delayMs ?? 600;
  await page.evaluate(async ({label, delayMs}) => {
    const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
    chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    await new Promise((r) => setTimeout(r, delayMs));
    const item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((m) => m.textContent!.trim() === label) as HTMLElement;
    (item.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, {label, delayMs});
}
