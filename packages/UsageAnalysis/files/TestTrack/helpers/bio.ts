import {Page} from '@playwright/test';

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
