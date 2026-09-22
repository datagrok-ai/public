/* The steps only the filter panel needs. Everything the panel itself reports — the cards, their
   criteria, the category and checkbox rows, the header counter — is the library's `viewers` tier
   over the `filter panel` element and the `filter card` kind (`grok-bdd list-steps`); what is left
   here is the counter's own menu, a card's search box and the hierarchical card, whose tree is not
   a grid and reports no hit areas: its rows are addressed by a "/"-separated path of the labels the
   product draws. */
import {expect, Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, exactText, gestures, viewers} from '@datagrok-libraries/bdd/runtime';

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

/** As a user pastes cells copied out of a table — the card reads a pasted list differently from
 * typed text. */
export const pasteIntoCardSearch = When('user pastes {string} into the search of the {string} filter card', async (page: Page, text: string, caption: string) => {
  await gestures.paste(page, cardSearch(page, caption), text);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'through the clipboard and the paste key into the card\'s search box; "\\n" is a line break'});

export const clearCardSearch = When('user clears the search of the {string} filter card', async (page: Page, caption: string) => {
  const input = cardSearch(page, caption);
  await input.click();
  await input.fill('');
  await viewers.settleAll(page);
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

// --- a second table on the same column object -------------------------------------------------------

/* GROK-13582 is about two tables holding one and the same column object: a table built on a copy of
   the column is independent by construction and could never show the leak. The platform has no
   gesture that builds such a table, so it is built through the API, and the sharing is proved — a
   tag written through the first table's column is read back through the second's — before the step
   ends. The new table gets a view of its own and becomes current. */
export const openSharedTable = Given('user opens a table {string} that shares the {string} column of the current table',
  async (page: Page, name: string, column: string) => {
    const shared = await page.evaluate(([n, c]) => {
      const w = window as any;
      const df = w.grok.shell.tv.dataFrame;
      const col = df.col(c);
      if (!col)
        throw new Error(`no "${c}" column in ${df.name}`);
      const other = w.DG.DataFrame.fromColumns([col]);
      other.name = n;
      w.grok.shell.addTableView(other);
      col.setTag('bdd-shared-probe', n);
      const back = other.col(c)?.getTag('bdd-shared-probe');
      col.setTag('bdd-shared-probe', null);
      return back === n;
    }, [name, column] as [string, string]);
    if (!shared)
      throw new Error(`"${name}" was built on a copy of "${column}", not on the same column object`);
    await page.waitForFunction((n) => (window as any).grok.shell.tv?.dataFrame?.name === n, name);
  }, {tier: 'api', description: 'a table built on the very column object of the current one (proved by a tag read back through it), with a view of its own'});

// --- links between tables ---------------------------------------------------------------------------

/* The graph a linked-tables scenario works on. Data > Link Tables builds a link one key column at a
   time through canvas column pickers, and the scenarios are about what travels over a link, not
   about building it, so the links are made through grok.data.linkTables; changing a link's type is
   done in the dialog, where the platform offers it. The link type is named as the dialog names it. */
export const linkTables = Given('the {string} table is linked to the {string} table by {string} to {string} as {string}',
  (page: Page, from: string, to: string, fromKeys: string, toKeys: string, type: string) => page.evaluate(([a, b, ka, kb, t]) => {
    const w = window as any;
    const table = (n: string): any => {
      const found = w.grok.shell.tables.find((x: any) => x.name === n);
      if (!found)
        throw new Error(`no table "${n}" is open; open: ${w.grok.shell.tables.map((x: any) => x.name).join(', ')}`);
      return found;
    };
    const sync = Object.values(w.DG.SYNC_TYPE).find((v) => v === t);
    if (!sync)
      throw new Error(`"${t}" is not a link type; the dialog offers: ${Object.values(w.DG.SYNC_TYPE).join(', ')}`);
    const keys = (s: string): string[] => s.split(/\s*,\s*/).filter((k) => k.length > 0);
    w.grok.data.linkTables(table(a), table(b), keys(ka), keys(kb), [sync]);
  }, [from, to, fromKeys, toKeys, type] as [string, string, string, string, string]),
  {tier: 'api', description: 'grok.data.linkTables with the key columns (comma-separated) of each side and the link type as the Link Tables dialog words it ("filter to filter")'});

// --- a short column picker in a dialog --------------------------------------------------------------

/* The Multi Value dialog's column selector offers only the table's categorical columns; a picker
   that short opens with no search box, so the library's typed pick has nothing to type into, and a
   click on a row of its list closes it without taking the row. The selector walks its columns on
   the arrow keys while it holds the focus, which is what a keyboard user does. */
export const pickColumnInDialog = When('user picks {string} in the column selector of the {string} dialog',
  async (page: Page, column: string, title: string) => {
    const dialog = page.locator('.d4-dialog').filter({has: page.locator('.d4-dialog-title', {hasText: exactText(title)})}).last();
    const selector = dialog.locator('.d4-column-selector').first();
    const shown = selector.locator('.d4-column-selector-column');
    await expect(selector, `the column selector of the "${title}" dialog`).toBeVisible({timeout: 5000});
    await selector.focus();
    const seen: string[] = [(await shown.textContent() ?? '').trim()];
    while (seen[seen.length - 1] !== column) {
      await selector.press('ArrowDown');
      const now = (await shown.textContent() ?? '').trim();
      if (now === seen[seen.length - 1] || seen.includes(now))
        throw new Error(`the column selector of the "${title}" dialog does not offer "${column}"; it offers: ${seen.join(', ')}`);
      seen.push(now);
    }
  }, {tier: 'ui', description: 'the selector takes the focus and the down arrow walks its columns until it shows this one; a column it never reaches fails, naming the ones it offered'});

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
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'the caption of one card dragged onto the top edge of another — the panel puts it before that card'});

// --- a histogram card's min and max fields --------------------------------------------------------

export const enterCardBound = When('user enters {string} into the {word} field of the {string} filter card',
  async (page: Page, value: string, which: string, caption: string) => {
    if (which !== 'min' && which !== 'max')
      throw new Error(`a histogram card has a min and a max field, not a "${which}" one`);
    const field = card(page, caption).locator(`input.d4-filter-input-${which}`).filter({visible: true}).first();
    await expect(field, `the ${which} field of the "${caption}" filter card (the card menu's "Min / max" shows it)`).toBeVisible({timeout: 5000});
    await field.click();
    await gestures.typeVerified(field, value, `the ${which} field of the "${caption}" filter card`);
    await field.press('Enter');
    await viewers.settleAll(page);
  }, {tier: 'ui', description: 'types the end into the field the card menu\'s "Min / max" reveals and commits it with Enter'});

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

/** The node box of a row: its hidden checkbox and the glyph the product draws in its place. */
const hierarchicalNode = (page: Page, path: string): Locator =>
  hierarchicalRow(page, path).locator('xpath=ancestor::div[contains(@class, "d4-tree-view-node")][1]');

export const clickHierarchicalCheckbox = When('user clicks on the checkbox of the {string} row of the hierarchical filter card', async (page: Page, path: string) => {
  const box = hierarchicalNode(page, path).locator('.d4-hierarchical-filter-checkbox-container').first();
  await box.click({timeout: 5000});
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'the box left of the row — it toggles that branch and leaves its siblings as they are, where a click on the label keeps the branch alone'});

/* The checkbox itself says nothing about a partial branch: the node's box carries `aria-checked`
   (true, false or mixed), set together with the glyph the tree draws. */
const ARIA_CHECKED: Record<string, string> = {checked: 'true', unchecked: 'false', partially: 'mixed'};

export const hierarchicalRowState = Then('the {string} row of the hierarchical filter card should be {word}( checked)',
  (page: Page, path: string, word: string) => {
    if (!ARIA_CHECKED[word])
      throw new Error(`a row of the hierarchical card is checked, unchecked or partially checked, not "${word}"`);
    return expect(hierarchicalNode(page, path).locator('.d4-hierarchical-filter-checkbox-container').first(),
      `the "${path}" row of the hierarchical filter card`).toHaveAttribute('aria-checked', ARIA_CHECKED[word], {timeout: 5000});
  }, {description: 'the aria-checked of the box left of the row: a branch whose children disagree says mixed'});

export const hierarchicalRowCounts = Then('the {string} row of the hierarchical filter card should count {int} rows',
  (page: Page, path: string, count: number) =>
    expect(hierarchicalRow(page, path).locator('.d4-hierarchical-filter-caption-count'),
      `the count of the "${path}" row of the hierarchical filter card`).toHaveText(String(count), {timeout: 5000}),
  {description: 'the number the row shows to its right — the rows of that branch that pass the filter'});
