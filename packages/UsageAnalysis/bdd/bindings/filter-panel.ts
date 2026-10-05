/* The steps only these features need. The filter panel's own vocabulary — the cards, their menus,
   search boxes, min/max fields and saved states, the counter's menu — is the library's `viewers`
   tier (`grok-bdd list-steps`); what is left here is a second table on the same column object and
   the hierarchical card, whose tree is not a grid and reports no hit areas: its rows are addressed
   by a "/"-separated path of the labels the product draws. */
import {Locator, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {exactText, expect, viewers} from '@datagrok-libraries/bdd/runtime';

const PANEL = '[name="viewer-Filters"]';

/** The filter panel of the view on screen (a cloned view leaves a second, hidden one behind). */
const panel = (page: Page): Locator => page.locator(PANEL).filter({visible: true}).first();

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
