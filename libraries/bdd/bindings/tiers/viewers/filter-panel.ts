import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from '../../../src/runtime/patience.js';
import {Given, Then, When} from '../../../src/registry.js';
import {el} from '../../../src/runtime/args.js';
import {atFeatureEnd} from '../../../src/runtime/harness.js';
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

/* Save or Apply > Save... keeps a named state on the server (a SavedFilter entity), where every later
   feature would find it: the user's own presets of that exact name are deleted by their ids. States
   kept in the browser before presets were entities ("filter-states") are dropped too, or the panel
   would move them to the server. */
const forgetPresets = (page: Page, names: string[]): Promise<void> => page.evaluate(async (all) => {
  const grok = (window as any).grok;
  const local = JSON.parse(localStorage.getItem('filter-states') ?? '{}');
  for (const n of all)
    delete local[n];
  localStorage.setItem('filter-states', JSON.stringify(local));
  const headers = {Authorization: localStorage.getItem('auth') ?? String(grok.dapi.token ?? '')};
  for (const n of all) {
    const found: {[key: string]: any}[] = await (await fetch(`${grok.dapi.root}/entities?text=${encodeURIComponent(n)}`, {headers})).json();
    for (const f of found.filter((e) => e['#type'] === 'SavedFilter' && e.friendlyName === n && e.author?.id === grok.shell.user.id))
      await fetch(`${grok.dapi.root}/domains/filters/${f.id}`, {method: 'DELETE', headers});
  }
}, names);

export const forgetSavedState = Given('no saved filter state {string} is kept, now or when the feature ends', async (page: Page, name: string) => {
  await forgetPresets(page, [name]);
  atFeatureEnd(page, () => forgetPresets(page, [name]));
}, {tier: 'api', description: 'the user\'s own saved filter presets of that name are deleted now, and again when the feature ends'});

/* Presets made through the server one after another, so the last is the latest: "#" in the name is
   its number. Each keeps a card on the text column, which every table with that column is offered. */
export const keepPresets = Given('{int} saved filter presets {string} on the text column {string} are kept until the feature ends',
  async (page: Page, count: number, pattern: string, column: string) => {
    const names = Array.from({length: count}, (_, i) => pattern.replace('#', String(i + 1)));
    await forgetPresets(page, names);
    atFeatureEnd(page, () => forgetPresets(page, names));
    await page.evaluate(async ([all, c]) => {
      const grok = (window as any).grok;
      const headers = {Authorization: localStorage.getItem('auth') ?? String(grok.dapi.token ?? ''), 'content-type': 'application/json'};
      for (const name of all) {
        const body = {name, states: {[c]: {type: 'categorical', column: c, columnType: 'string'}}, columns: [`${c}:string`]};
        const answer = await fetch(`${grok.dapi.root}/domains/filters`, {method: 'POST', headers, body: JSON.stringify(body)});
        if (!answer.ok)
          throw new Error(`saving the preset "${name}": ${answer.status} ${await answer.text()}`);
      }
    }, [names, column] as [string[], string]);
  }, {tier: 'api', description: 'saved one after another, so the last is the latest; they are deleted by name now and when the feature ends'});

/* A preset's row under Save or Apply (opened first by the steps that act on it): its name, "who, when" beside it, and while hovered its share and
   delete icons in that place. */
const presetRow = (page: Page, name: string): Locator => page.locator('.d4-menu-popup .d4-filter-preset').filter({visible: true})
  .filter({has: page.locator(':scope > .d4-menu-item-label', {hasText: exactText(name)})}).first();
const ICONS: {[name: string]: string} = {share: 'fa-share-alt', delete: 'fa-trash-alt'};

export const presetBy = Then('the filter preset {string} in the open menu should be by {string}', async (page: Page, name: string, who: string) => {
  await expect.poll(async () => (await presetRow(page, name).locator(':scope > .d4-menu-item-shortcut').textContent().catch(() => null))?.split(',')[0],
    {message: `who saved the "${name}" preset, beside its name`}).toBe(who);
}, {description: 'the author beside its name ("me" for your own), before when it was last saved'});

export const presetOffers = Then('the filter preset {string} in the open menu should offer {string}', async (page: Page, name: string, list: string) => {
  await v.menuLabels(page, 'Save or Apply');
  const row = presetRow(page, name);
  await row.hover();
  const wanted = list === 'nothing' ? [] : list.split(/\s*(?:,|and)\s*/);
  await expect.poll(() => row.locator('.d4-filter-preset-actions .grok-icon').filter({visible: true})
    .evaluateAll((icons, all) => Object.keys(all).filter((k) => icons.some((i) => i.classList.contains(all[k]))), ICONS),
  {message: `the icons of the hovered "${name}" preset`}).toEqual(wanted);
}, {tier: 'ui', description: '"share and delete", "share" or "nothing": the icons that show while the row is hovered'});

export const clickPresetIcon = When('user clicks the {word} icon of the filter preset {string} in the open menu', async (page: Page, icon: string, name: string) => {
  if (!ICONS[icon])
    throw new Error(`a preset row offers ${Object.keys(ICONS).join(' and ')}, not "${icon}"`);
  await v.menuLabels(page, 'Save or Apply');
  const row = presetRow(page, name);
  await row.hover();
  await row.locator(`.d4-filter-preset-actions .${ICONS[icon]}`).click();
}, {tier: 'ui', description: 'share opens the Share dialog, delete asks to confirm; either closes the menu'});

export const presetNamesAuthor = Then('the filter preset {string} in the open menu should name its author', async (page: Page, name: string) => {
  await v.menuLabels(page, 'Save or Apply');
  await expect.poll(async () => (await presetRow(page, name).locator(':scope > .d4-menu-item-shortcut').textContent().catch(() => null))?.split(',')[0].trim(),
    {message: `who saved the "${name}" preset, beside its name`}).toMatch(/^(?!me$).+/);
}, {description: 'someone else\'s preset shows its author\'s name beside it, where one\'s own shows "me"'});

export const scrollPresets = When('user scrolls the presets under Save or Apply to the end in the open menu', async (page: Page) => {
  await v.menuLabels(page, 'Save or Apply');
  const list = page.locator('.d4-filter-presets-menu').filter({visible: true});
  const box = await list.boundingBox();
  if (!box)
    throw new Error('the Save or Apply menu is not open');
  await page.mouse.move(box.x + box.width / 2, box.y + box.height - 20);
  const before = await list.evaluate((e) => e.scrollHeight);
  await page.mouse.wheel(0, before);
  await expect.poll(() => list.evaluate((e) => e.scrollHeight), {message: 'the next presets added below'}).toBeGreaterThan(before);
}, {tier: 'ui', description: 'the wheel over the preset list takes it to its end, which adds the next 20 presets below'});

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

// the header's column picker clears its own label, so the base `select` cannot wait for one
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

// opened on blank space: a right-click in the middle of a card belongs to that card's grid
export const pickPanelMenu = When('user picks {string} from the filter panel menu', async (page: Page, path: string) => {
  const areas = await v.hitAreas(page, el('filter panel'));
  const view = areas['view'];
  if (!view)
    throw new Error('the filter panel reports no "view" area');
  const bottom = Object.keys(areas).filter((k) => k.startsWith('card '))
    .reduce((y, k) => Math.max(y, areas[k].y + areas[k].height), view.y);
  if (bottom <= view.y + view.height - 8)
    await v.openContextMenuAt(page, view.x + view.width / 2, (bottom + view.y + view.height) / 2);
  else {
    const dock = panel(page).locator('xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]');
    const icon = dock.locator('.panel-titlebar [name="icon-font-icon-menu"]').first();
    if (await icon.count() === 0)
      throw new Error('the cards fill the filter panel and its title bar shows no menu icon (not docked, or the dockShowMenuLink setting off): nowhere to open the panel\'s own menu from');
    await g.hover(page, el('filter panel'));
    await icon.click();
  }
  await v.pickMenuPath(page, path);
}, {tier: 'ui', description: 'a right-click on the panel below its last card — the cards own the rest of it — or, when the cards fill it, the menu icon of its docked title bar (the same menu); then the menu path ("Add Filter | Hierarchical")'});

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
