/* The `viewers` tier, the filter panel's entry points that have no element to click: the header's
   column picker clears its own label, so the base `select` cannot wait for one, and the panel's own
   context menu has to be opened on blank space, since a right-click in the middle of a card belongs
   to that card's grid. Everything the panel reports — cards, criteria, the header counter — is the
   `filter panel` element and the `filter card` kind. */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from '../../../src/runtime/patience.js';
import {Then, When} from '../../../src/registry.js';
import {el} from '../../../src/runtime/args.js';
import {atFeatureEnd} from '../../../src/runtime/harness.js';
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

/* --- every view's filter panel, across a project round trip ------------------------------------------
   Per table view, the cards that filter and their summaries, read from each panel's own readings. A
   Scaffold Tree card (Chem) has no caption among the panel's readings: its tree says what it filters
   by. Not the table's own list of filters (`df.rows.filters`): a panel of a view restored with a
   project joins it only once it is attached, so it holds what was shown, not what was saved. Kept in
   the test run, not in the page, and forgotten when the feature ends. */
async function panelStates(page: Page): Promise<Record<string, string>> {
  return page.evaluate(() => {
    const w = window as any;
    const out: Record<string, string> = {};
    for (const tv of Array.from(w.grok.shell.tableViews as any[])) {
      const host = Array.from((tv.root as HTMLElement).querySelectorAll('[name="viewer-Filters"]'))[0];
      if (!host) {
        out[tv.name] = 'no panel';
        continue;
      }
      const viewer = w.__bdd?.viewerOf(host) ?? w.DG.Widget.find(host);
      const values = viewer?.getWidgetStatus?.()?.values ?? {};
      const cards = String(values['cards'] ?? '').split(', ').filter((c) => values[`filtering of ${c}`] === true);
      const trees = (tv.getFiltersGroup().filters as any[]).filter((f) => f?.viewer && 'checkedNodes' in f.viewer)
        .map((f) => {
          const t = f.viewer.getWidgetStatus().values;
          return `scaffold tree [${t['checked nodes']} of ${t['nodes']} checked, ${t['bit operation']}]`;
        });
      out[tv.name] = `${values['active'] === false ? 'off: ' : ''}${[...cards.map((c) => `${c} [${values[`summary of ${c}`] ?? ''}]`), ...trees].join('; ')}`;
    }
    return out;
  });
}

const rememberedPanels = new WeakMap<Page, Record<string, string>>();

export const rememberPanelStates = When('user remembers what the filter panel of every view filters by', async (page: Page) => {
  await v.settleAll(page);
  if (!rememberedPanels.has(page))
    atFeatureEnd(page, async () => { rememberedPanels.delete(page); });
  rememberedPanels.set(page, await panelStates(page));
}, {tier: 'api', description: 'per table view, the cards that filter and their summaries — the state a saved project has to bring back'});

export const panelStatesAsRemembered = Then('the filter panel of every view should filter by what was remembered', async (page: Page) => {
  const want = rememberedPanels.get(page);
  if (!want)
    throw new Error('nothing was remembered — "user remembers what the filter panel of every view filters by" comes first');
  await expect.poll(() => panelStates(page), {message: 'what the filter panel of each view filters by'}).toEqual(want);
}, {description: 'every table view has the same filtering cards with the same summaries as remembered'});

/** A substructure card searches in the background ("searching of <col>"), and a Scaffold Tree counts
 * the hits of its nodes after it loads (-1 until then): a count read before both are done is the count
 * of a filter not applied yet. A substructure card also lays itself out after the panel lists it
 * ("drawing of <col>", reported once it is on the page), and the cards below it move until it has —
 * in the current view: a view restored with a project and not shown since has cards not laid out yet. */
export const filtersDoneComputing = Then('the filters of every view should have finished computing', async (page: Page) => {
  let pending: string[] = [];
  await expect.poll(async () => (pending = await page.evaluate(() => {
    const w = window as any;
    const out: string[] = [];
    for (const tv of Array.from(w.grok.shell.tableViews as any[])) {
      const host = (tv.root as HTMLElement).querySelector('[name="viewer-Filters"]');
      if (!host)
        continue;
      const values = w.__bdd?.viewerOf(host)?.getWidgetStatus?.()?.values ?? {};
      for (const [k, s] of Object.entries(values))
        if (k.startsWith('searching of ') && s === true)
          out.push(`${tv.name}: ${k}`);
      const entries = Object.entries(values);
      const structureCards = entries.filter(([k, s]) => k.startsWith('type of ') && s === 'Chem:substructureFilter').length;
      const drawn = entries.filter(([k, s]) => k.startsWith('drawing of ') && s === false).length;
      if (tv === w.grok.shell.tv && drawn < structureCards)
        out.push(`${tv.name}: ${structureCards - drawn} structure card(s) drawing`);
      for (const f of tv.getFiltersGroup().filters as any[]) {
        if (!(f?.viewer && 'checkedNodes' in f.viewer))
          continue;
        const tree = f.viewer.getWidgetStatus().values;
        for (let i = 1; i <= Number(tree['nodes'] ?? 0); i++)
          if (tree[`checked of node ${i}`] === true && Number(tree[`hits of node ${i}`]) < 0)
            out.push(`${tv.name}: hits of checked node ${i}`);
      }
    }
    return out;
  })).length, {timeout: pollMs(60000), message: 'filters still computing'}).toBe(0).catch(() => {
    throw new Error(`filters still computing: ${pending.join(', ')}`);
  });
  await v.settleAll(page);
}, {description: 'no substructure card of any view is searching, none of the current view is still drawing, and every checked Scaffold Tree node has counted its hits (up to a minute)'});
