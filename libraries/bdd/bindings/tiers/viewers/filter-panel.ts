/* The `viewers` tier, the filter panel's entry points that have no element to click: the header's
   column picker clears its own label, so the base `select` cannot wait for one, and the panel's own
   context menu has to be opened on blank space, since a right-click in the middle of a card belongs
   to that card's grid. Everything the panel reports — cards, criteria, the header counter — is the
   `filter panel` element and the `filter card` kind. */
import {Locator, Page} from '@playwright/test';
import {expect} from '../../../src/runtime/patience.js';
import {Given, Then, When} from '../../../src/registry.js';
import {atFeatureEnd} from '../../../src/runtime/harness.js';
import {el} from '../../../src/runtime/args.js';
import {exactText} from '../../../src/runtime/locate.js';
import * as g from '../../../src/runtime/gestures.js';
import * as v from '../../../src/runtime/viewers.js';

/** The filter panel of the view on screen (a cloned view leaves a second, hidden one behind). */
const panel = (page: Page): Locator => page.locator('[name="viewer-Filters"]').filter({visible: true}).first();

// --- the counter's own menu ------------------------------------------------------------------

export const pickCounterMenu = When('user picks {string} from the filter counter menu', async (page: Page, path: string) => {
  await g.hover(page, el('filter panel'));
  await panel(page).locator('[name="active-filter-counter"]').first().click();
  await v.pickMenuPath(page, path);
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
    await v.pickMenuPath(page, path);
  }, {tier: 'ui', description: 'a path in the card\'s own menu ("Mode | Radio")'});

export const typeIntoCardSearch = When('user types {string} into the search of the {string} filter card', async (page: Page, text: string, caption: string) => {
  const input = cardSearch(page, caption);
  await input.click();
  await input.fill(text);
  await v.settleAll(page);
}, {tier: 'ui', description: 'the search box the card\'s search icon reveals — the categorical card narrows its rows to the matches, the hierarchical one hides the nodes that do not match'});

/** As a user pastes cells copied out of a table — the card reads a pasted list differently from
 * typed text. */
export const pasteIntoCardSearch = When('user pastes {string} into the search of the {string} filter card', async (page: Page, text: string, caption: string) => {
  await g.paste(page, cardSearch(page, caption), text);
  await v.settleAll(page);
}, {tier: 'ui', description: 'through the clipboard and the paste key into the card\'s search box; "\\n" is a line break'});

export const clearCardSearch = When('user clears the search of the {string} filter card', async (page: Page, caption: string) => {
  const input = cardSearch(page, caption);
  await input.click();
  await input.fill('');
  await v.settleAll(page);
}, {tier: 'ui'});

// --- the panel's saved states ------------------------------------------------------------------------

/* Save or Apply > Save... keeps a named state in the browser's localStorage under "filter-states"
   (FiltersCore.FILTER_STATES_LS_KEY), where every later feature on the worker's page would find it. */
export const forgetSavedState = Given('no saved filter state {string} is kept, now or when the feature ends', async (page: Page, name: string) => {
  const forget = (): Promise<void> => page.evaluate((n) => {
    const all = JSON.parse(localStorage.getItem('filter-states') ?? '{}');
    delete all[n];
    localStorage.setItem('filter-states', JSON.stringify(all));
  }, name);
  await forget();
  atFeatureEnd(page, forget);
}, {tier: 'api', description: 'the state of that name is taken out of the browser\'s saved filter states now, and again when the feature ends'});

// --- moving a card --------------------------------------------------------------------------------

/* A card is moved by its caption: the press starts the platform's drag, the pointer walks in small
   steps (a single jump never passes the start threshold) to the top edge of the other card —
   the strip the panel drops into — and the release there puts it before that card. */
export const dragCard = When('user drags the {string} filter card above the {string} filter card',
  async (page: Page, moved: string, other: string) => {
    const from = await card(page, moved).locator('.d4-filter-column-name').boundingBox();
    const target = await card(page, other).boundingBox();
    if (!from || !target)
      throw new Error(`the "${from ? other : moved}" filter card is not on screen`);
    const a = {x: from.x + Math.min(10, from.width / 2), y: from.y + from.height / 2};
    const b = {x: a.x, y: target.y + 3};
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(b.x, b.y, {steps: 12});
    await page.mouse.up();
    await v.settleAll(page);
  }, {tier: 'ui', description: 'the caption of one card dragged onto the top edge of another — the panel puts it before that card'});

// --- a histogram card's min and max fields --------------------------------------------------------

export const enterCardBound = When('user enters {string} into the {word} field of the {string} filter card',
  async (page: Page, value: string, which: string, caption: string) => {
    if (which !== 'min' && which !== 'max')
      throw new Error(`a histogram card has a min and a max field, not a "${which}" one`);
    const field = card(page, caption).locator(`input.d4-filter-input-${which}`).filter({visible: true}).first();
    await expect(field, `the ${which} field of the "${caption}" filter card (the card menu's "Min / max" shows it)`).toBeVisible({timeout: 5000});
    await field.click();
    await g.typeVerified(field, value, `the ${which} field of the "${caption}" filter card`);
    await field.press('Enter');
    await v.settleAll(page);
  }, {tier: 'ui', description: 'types the end into the field the card menu\'s "Min / max" reveals and commits it with Enter'});


// --- the substructure card (Chem) -----------------------------------------------------------------

/* The substructure filter card's own controls: the search-type choice under the sketch area and, once
   the card's settings icon is on, the fingerprint choice and the similarity cutoff. They carry no
   label (the search type) or a label the platform does not name (FP, the slider), so no input kind
   reaches them. What the card holds is read on the filter panel: "search type of <column>" and the
   rest of the card's readings. */
const searchType = (page: Page, caption: string): Locator => card(page, caption).locator('.chem-filter-search-type select').first();

export const pickSearchType = When('user picks search type {string} in the {string} filter card', async (page: Page, type: string, caption: string) => {
  const select = searchType(page, caption);
  await select.waitFor({state: 'visible', timeout: 10000});
  await select.selectOption({label: type});
  await v.settleAll(page);
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
  await v.settleAll(page);
}, {tier: 'ui', description: 'the number box beside the cutoff slider, shown for the Similar search type'});

export const addCardFor = When('user adds a card for {string} to the filter panel', async (page: Page, column: string) => {
  await g.hover(page, el('filter panel'));
  const selector = panel(page).locator('[name="div-column-combobox-add-filter"]').first();
  const box = await selector.boundingBox();
  if (!box)
    throw new Error('the filter panel shows no add-filter selector: its header controls appear on hover, and the panel is not hovered');
  await page.mouse.move(box.x + Math.min(10, box.width / 2), box.y + box.height / 2);
  await page.mouse.down();
  await page.mouse.up();
  // the plus icon restores whatever was focused before the click, so the selector — which is what
  // the typed name goes to — loses the focus its own mouse-down gave it. The pointer stays on the
  // panel: the header this picker belongs to is only shown while the panel is hovered
  await g.pickInColumnGrid(page, column, 'the filter panel', selector);
  await panel(page).locator('.d4-filter')
    .filter({has: page.locator('.d4-filter-column-name', {hasText: exactText(column)})})
    .first().waitFor({state: 'visible'});
  await v.settleAll(page);
}, {tier: 'ui', description: 'the plus selector in the panel header: opens its column grid, types the name and commits — the card appears at the top of the panel'});

export const pickPanelMenu = When('user picks {string} from the filter panel menu', async (page: Page, path: string) => {
  const areas = await v.hitAreas(page, el('filter panel'));
  const view = areas['view'];
  if (!view)
    throw new Error('the filter panel reports no "view" area');
  const bottom = Object.keys(areas).filter((k) => k.startsWith('card '))
    .reduce((y, k) => Math.max(y, areas[k].y + areas[k].height), view.y);
  if (bottom > view.y + view.height - 8)
    throw new Error('the cards fill the filter panel: no blank space left to open the panel\'s own menu on');
  await v.openContextMenuAt(page, view.x + view.width / 2, (bottom + view.y + view.height) / 2);
  await v.pickMenuPath(page, path);
}, {tier: 'ui', description: 'a right-click on the panel below its last card — the cards own the rest of it — then the menu path ("Add Filter | Hierarchical")'});

export const panelHasNoCardOfType = Then('the filter panel should have no {string} filter card', async (page: Page, type: string) => {
  let cards: string[] = [];
  await expect.poll(async () => {
    const values: Record<string, unknown> = await v.onViewer(page, el('filter panel'), (e) => (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {});
    cards = Object.entries(values).filter(([k, t]) => k.startsWith('type of ') && t === type).map(([k]) => k.slice(8));
    return cards.length;
  }, {message: `${type} cards on the filter panel`}).toBe(0).catch(() => {
    throw new Error(`the filter panel still has ${type} cards on: ${cards.join(', ')}`);
  });
}, {description: 'no card whose "type of <caption>" reading is the given filter type ("Chem:substructureFilter")'});
