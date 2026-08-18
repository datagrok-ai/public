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
    (item.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, {label, delayMs});
}
