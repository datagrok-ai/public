/* The substructure filter card's own controls: the search-type choice under the sketch area and, once
   the card's settings icon is on, the fingerprint choice and the similarity cutoff. They carry no
   label (the search type) or a label the platform does not name (FP, the slider), so no input kind
   reaches them. What the card holds is read on the filter panel: "search type of <column>" and the
   rest of the card's readings. */
import {Locator, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {exactText, expect, viewers} from '@datagrok-libraries/bdd/runtime';

const card = (page: Page, caption: string): Locator => page.locator('[name="viewer-Filters"]').filter({visible: true}).first()
  .locator('.d4-filter').filter({has: page.locator('.d4-filter-column-name', {hasText: exactText(caption)})}).first();

const searchType = (page: Page, caption: string): Locator => card(page, caption).locator('.chem-filter-search-type select').first();

export const pickSearchType = When('user picks search type {string} in the {string} filter card', async (page: Page, type: string, caption: string) => {
  const select = searchType(page, caption);
  await select.waitFor({state: 'visible', timeout: 10000});
  await select.selectOption({label: type});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the choice under the card\'s sketch area ("Contains", "Not contains", "Exact", ...)'});

export const offersSearchTypes = Then('the {string} filter card should offer search types {string}', async (page: Page, caption: string, list: string) => {
  const select = searchType(page, caption);
  await expect.poll(async () => (await select.locator('option').allTextContents()).map((o) => o.trim()),
    {message: `the search types of the "${caption}" card`}).toEqual(list.split(/\s*,\s*/));
}, {description: 'the options of the search-type choice, in order'});

export const openCardSettings = When('user opens the settings of the {string} filter card', async (page: Page, caption: string) => {
  const select = searchType(page, caption);
  if (await select.isVisible())
    return;
  await card(page, caption).hover();
  await card(page, caption).locator('.chem-search-options-icon').first().click();
  await select.waitFor({state: 'visible', timeout: 5000});
}, {tier: 'ui', description: 'the gear icon of the card, which toggles the search type, fingerprint and cutoff controls; left alone when they already show'});

export const setCutoff = When('user sets the similarity cutoff of the {string} filter card to {float}', async (page: Page, caption: string, value: number) => {
  const editor = card(page, caption).locator('.chem-filter-sim-cutoff-editor').first();
  await editor.waitFor({state: 'visible', timeout: 5000});
  await editor.fill(String(value));
  await editor.press('Enter');
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the number box beside the cutoff slider, shown for the Similar search type'});
