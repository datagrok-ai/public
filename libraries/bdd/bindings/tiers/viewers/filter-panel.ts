/* The `viewers` tier, the filter panel's entry points that have no element to click: the header's
   column picker clears its own label, so the base `select` cannot wait for one, and the panel's own
   context menu has to be opened on blank space, since a right-click in the middle of a card belongs
   to that card's grid. Everything the panel reports — cards, criteria, the header counter — is the
   `filter panel` element and the `filter card` kind. */
import {Locator, Page} from '@playwright/test';
import {expect} from '../../../src/runtime/patience.js';
import {Then, When} from '../../../src/registry.js';
import {el} from '../../../src/runtime/args.js';
import {exactText} from '../../../src/runtime/locate.js';
import * as g from '../../../src/runtime/gestures.js';
import * as v from '../../../src/runtime/viewers.js';

/** The filter panel of the view on screen (a cloned view leaves a second, hidden one behind). */
const panel = (page: Page): Locator => page.locator('[name="viewer-Filters"]').filter({visible: true}).first();

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
