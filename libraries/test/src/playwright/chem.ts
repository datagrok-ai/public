/**
 * Playwright helpers shared by Chem/* spec files.
 *
 * Each helper is a verbatim extraction of a block previously pasted into the
 * Chem specs — same selectors, same mechanism, same sleeps. Imported as:
 * `import * as chem from '@datagrok-libraries/test/src/playwright/chem';` (or `'@datagrok-libraries/test/src/playwright/chem'` from
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
    const menuItem = item.closest('.d4-menu-item') as HTMLElement;
    // One click on [name="div-Chem"] puts every menubar root's labels in the document —
    // 165 of them on dev — so an exact text match can land on another feature's item and
    // click it silently. Three leaves collide repo-wide: Activity Cliffs..., Hierarchical
    // Clustering... (ML's copy comes first in DOM order), Boltz... For those, address the
    // item by its full-path name instead: [name="div-Chem---Analyze---Hierarchical-Clustering..."].
    const named = item.closest('[name^="div-"]');
    const owner = named ? named.getAttribute('name') : null;
    if (owner !== null && !owner.startsWith('div-Chem'))
      throw new Error(`openChemMenuItem("${label}") matched ${owner}, which is not a Chem menu item`);
    menuItem.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, {label, delayMs});
}
