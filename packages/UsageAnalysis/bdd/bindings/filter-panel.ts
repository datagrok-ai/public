/* The steps only the filter panel needs. Everything the panel itself reports — the cards, their
   criteria, the category and checkbox rows, the header counter — is the library's `viewers` tier
   over the `filter panel` element and the `filter card` kind (`grok-bdd list-steps`); what is left
   here is the counter's own menu, a card's search box and the hierarchical card, whose tree is not
   a grid and reports no hit areas: its rows are addressed by a "/"-separated path of the labels the
   product draws. */
import {expect, Locator, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, exactText, gestures, viewers} from '@datagrok-libraries/bdd/runtime';

const PANEL = '[name="viewer-Filters"]';

/** The filter panel of the view on screen (a cloned view leaves a second, hidden one behind). */
const panel = (page: Page): Locator => page.locator(PANEL).filter({visible: true}).first();

// --- the counter's own menu ------------------------------------------------------------------

export const pickCounterMenu = When('user picks {string} from the filter counter menu', async (page: Page, path: string) => {
  await gestures.hover(page, el('filter panel'));
  await panel(page).locator('[name="active-filter-counter"]').first().click();
  await viewers.pickMenuPath(page, path);
}, {tier: 'ui', description: 'the counter opens its own menu on a left click — its tooltip says so, and a right-click there brings up the viewer\'s menu instead'});

// --- a card's own search box ---------------------------------------------------------------------

/** A card of the panel by the caption it shows. */
const card = (page: Page, caption: string): Locator => panel(page).locator('.d4-filter')
  .filter({has: page.locator('.d4-filter-column-name', {hasText: exactText(caption)})}).first();

/** The card's search box: a bare input the product appends to the card body, so no input kind
 * reaches it, and the card's `search` part finds the "Filter by search results" checkbox first. */
const cardSearch = (page: Page, caption: string): Locator => card(page, caption).locator('input.d4-search-input').first();

/** The card's own menu button: hover-revealed, and it opens its menu on a plain click. */
async function openIndicatorMenu(page: Page, caption: string): Promise<void> {
  await card(page, caption).hover();
  await card(page, caption).locator('.d4-filter-indicator').first().click();
  // without this a menu that never opened reads as a menu with nothing in it, and every
  // "should not contain" claim about it passes
  await page.locator('.d4-menu-popup').filter({visible: true}).first().waitFor({state: 'visible'});
}

export const openCardIndicatorMenu = When('user opens the indicator menu of the {string} filter card', (page: Page, caption: string) =>
  openIndicatorMenu(page, caption), {tier: 'ui', description: 'the counter at the right of the card caption — it shows while the card is hovered and opens the card\'s batch and mode menu'});

export const pickCardIndicatorMenu = When('user picks {string} from the indicator menu of the {string} filter card',
  async (page: Page, path: string, caption: string) => {
    await openIndicatorMenu(page, caption);
    await viewers.pickMenuPath(page, path);
  }, {tier: 'ui', description: 'a path in the card\'s own menu ("Mode | Radio")'});

export const typeIntoCardSearch = When('user types {string} into the search of the {string} filter card', async (page: Page, text: string, caption: string) => {
  const input = cardSearch(page, caption);
  await input.click();
  await input.fill(text);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the search box the card\'s search icon reveals — the categorical card narrows its rows to the matches, the hierarchical one hides the nodes that do not match'});

export const clearCardSearch = When('user clears the search of the {string} filter card', async (page: Page, caption: string) => {
  const input = cardSearch(page, caption);
  await input.click();
  await input.fill('');
  await viewers.settleAll(page);
}, {tier: 'ui'});

// --- the hierarchical card ---------------------------------------------------------------------

/** The card whose body is a hierarchy tree. */
const hierarchicalCard = (page: Page): Locator => panel(page).locator('.d4-filter')
  .filter({has: page.locator('.d4-hierarchical-filter-caption-value')}).first();

/** The label of the row a "/"-separated path names: every segment is looked up among the rows the
 * previous one opened, so "F / Caucasian" is the Caucasian under F and not the one under M. */
function hierarchicalRow(page: Page, path: string): Locator {
  const segments = path.split(/\s*\/\s*/).filter((s) => s.length > 0);
  let scope: Locator = hierarchicalCard(page);
  let label = scope.locator('.d4-hierarchical-filter-caption-label').first();
  for (let i = 0; i < segments.length; i++) {
    label = scope.locator('.d4-hierarchical-filter-caption-label')
      .filter({has: page.locator('.d4-hierarchical-filter-caption-value', {hasText: exactText(segments[i])})}).first();
    if (i < segments.length - 1)
      scope = label.locator('xpath=ancestor::div[contains(@class, "d4-tree-view-group")][1]')
        .locator('.d4-tree-view-group-host').first();
  }
  return label;
}

export const clickHierarchicalRow = When('user clicks on the {string} row of the hierarchical filter card', async (page: Page, path: string) => {
  await hierarchicalRow(page, path).click({timeout: 5000});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the row label, which the product reads as "this branch alone": it clears every other node and checks this one'});

export const expandHierarchicalRow = When('user expands the {string} row of the hierarchical filter card', async (page: Page, path: string) => {
  await hierarchicalRow(page, path).locator('xpath=ancestor::div[contains(@class, "d4-tree-view-group")][1]')
    .locator('.d4-tree-view-tri').first().click({timeout: 5000});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the twistie of a branch — the children of a level are built when it is first opened'});

export const hierarchicalListsRow = Then('the hierarchical filter card should list the {string} row', (page: Page, path: string) =>
  expect(hierarchicalRow(page, path), `the "${path}" row of the hierarchical filter card`).toBeVisible({timeout: 5000}));

// the tree keeps the nodes it has built and hides them, so "not listed" is hidden or never built
export const hierarchicalHidesRow = Then('the hierarchical filter card should not list the {string} row', (page: Page, path: string) =>
  expect(hierarchicalRow(page, path), `the "${path}" row of the hierarchical filter card`).toBeHidden({timeout: 5000}));

export const hierarchicalRowCounts = Then('the {string} row of the hierarchical filter card should count {int} rows',
  (page: Page, path: string, count: number) =>
    expect(hierarchicalRow(page, path).locator('.d4-hierarchical-filter-caption-count'),
      `the count of the "${path}" row of the hierarchical filter card`).toHaveText(String(count), {timeout: 5000}),
  {description: 'the number the row shows to its right — the rows of that branch that pass the filter'});
