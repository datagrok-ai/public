import {Page} from '@playwright/test';

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
