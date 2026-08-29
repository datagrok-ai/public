/* The control inspector on the live demo: clicking a u2 control inside U2Demo's content host makes
   it the current object (capture phase, non-freezing), and the u2-control ObjectHandler renders a
   PropertyGrid in the context panel. The chip check is the regression guard for that ordering: an
   entity chip's own bubble-phase `grok.shell.o = entity` freezes the channel and must win over the
   inspector's earlier write, so the panel shows the entity, not the chip. */
import {ok, shot} from '../local.mjs';
import {openDemoPage} from '../lib.mjs';

const CONTENT = '.u2demo-content';
const TEXT_INPUT = `${CONTENT} [data-u2="text-input"]`;
const GRID = '.u2-inspector [data-u2="property-grid"]';

export async function fixture(page) {
  await openDemoPage(page, 'basic-inputs');
  await page.waitForSelector(TEXT_INPUT, {timeout: 15000});
}

const gridShown = (page, timeout) =>
  page.waitForSelector(GRID, {timeout}).then(() => true).catch(() => false);

async function checkControlClick(page) {
  // the platform panel renders one `grok.shell.o` write behind at times — re-assert with the
  // click a user would make (same discipline as lib.mjs selectRow)
  await page.locator(TEXT_INPUT).first().click();
  let shown = await gridShown(page, 6000);
  if (!shown) {
    await page.locator(TEXT_INPUT).first().click();
    shown = await gridShown(page, 6000);
  }
  await shot(page, 'u2demo-inspector-1-grid');
  ok('u2demo-inspector/1a/clicking-a-text-input-renders-the-property-grid-in-the-context-panel',
    shown, `grid=${shown}`);

  const header = await page.evaluate(() =>
    document.querySelector('.u2-inspector h3')?.textContent ?? '');
  ok('u2demo-inspector/1b/the-panel-header-names-the-clicked-control-tag',
    header === 'u2-text-input', `header="${header}"`);
}

async function checkChipShowsEntity(page) {
  await page.locator('.u2demo-nav .u2-tree-row')
    .filter({has: page.locator('.u2-tree-label', {hasText: 'Entities'})}).first().click();
  await page.waitForSelector(`${CONTENT} .u2-chip`, {timeout: 15000});
  await page.locator(`${CONTENT} .u2-chip`).first().click();
  const settled = await page.waitForFunction((grid) => {
    const o = grok.shell.o;
    return o != null && o.dart != null && !(o.root instanceof HTMLElement) &&
      document.querySelector(grid) == null;
  }, GRID, {timeout: 10000}).then(() => true).catch(() => false);
  const detail = await page.evaluate(() => {
    const o = grok.shell.o;
    return {ctor: o?.constructor?.name ?? 'null', dart: o?.dart != null,
      inspectorStillShown: document.querySelector('.u2-inspector') != null};
  });
  await shot(page, 'u2demo-inspector-2-chip');
  ok('u2demo-inspector/2a/an-entity-chip-click-puts-the-entity-not-the-chip-control-in-the-panel',
    settled, `o=${detail.ctor} dart=${detail.dart} inspectorStillShown=${detail.inspectorStillShown}`);
}

/** What the inspector panel says, whatever state it is in. */
const panelSays = (page, text) => page.waitForFunction((t) =>
  (document.querySelector('.u2-inspector')?.textContent ?? '').includes(t), text, {timeout: 8000})
  .then(() => true).catch(() => false);

/** A click that hits no u2 control is a no-op, not a wipe: the panel keeps whatever it holds. The
 * unbounded walk this replaced fell through to the app's own splitter and rendered 15 blank rows;
 * the empty state that replaced THAT ate the source panel every time a heading was clicked. */
async function checkWhitespace(page) {
  // the chip check left the Entities leaf open; the panel has to hold a control grid to prove a
  // whitespace click does not take it away
  await page.locator('.u2demo-nav .u2-tree-row')
    .filter({has: page.locator('.u2-tree-label', {hasText: 'Basic inputs'})}).first().click();
  await page.waitForSelector(TEXT_INPUT, {timeout: 15000});
  await page.locator(TEXT_INPUT).first().click();
  let before = await gridShown(page, 6000);
  if (!before) {
    await page.locator(TEXT_INPUT).first().click();
    before = await gridShown(page, 6000);
  }
  // the content host's own padding: no control under the pointer
  await page.locator(CONTENT).click({position: {x: 3, y: 3}});
  await page.waitForTimeout(1200);
  const after = await page.evaluate((grid) => ({
    grid: document.querySelector(grid) != null,
    emptied: (document.querySelector('.u2-inspector')?.textContent ?? '')
      .includes('not a u2 component'),
  }), GRID);
  await shot(page, 'u2demo-inspector-3-whitespace');
  ok('u2demo-inspector/3a/clicking-bare-whitespace-leaves-the-panel-alone-instead-of-emptying-it',
    before && after.grid && !after.emptied,
    `gridBefore=${before} gridAfter=${after.grid} emptyState=${after.emptied}`);
}

async function checkUnregisteredControl(page) {
  await page.locator('.u2demo-nav .u2-tree-row')
    .filter({has: page.locator('.u2-tree-label', {hasText: 'Trees'})}).first().click();
  await page.waitForSelector(`${CONTENT} [data-u2="tree"]`, {timeout: 15000});
  await page.locator(`${CONTENT} [data-u2="tree"] .u2-tree-label`).first().click();
  let shown = await panelSays(page, 'No inspectable properties');
  if (!shown) {
    await page.locator(`${CONTENT} [data-u2="tree"] .u2-tree-label`).first().click();
    shown = await panelSays(page, 'No inspectable properties');
  }
  const text = await page.evaluate(() => document.querySelector('.u2-inspector')?.textContent ?? '');
  await shot(page, 'u2demo-inspector-4-unregistered');
  ok('u2demo-inspector/4a/a-control-the-registry-does-not-describe-says-so-instead-of-an-empty-grid',
    shown, `panel="${text.slice(0, 120)}"`);
}

export const checks = [
  {id: 'u2demo-inspector/1 a control click renders its properties in the context panel', run: checkControlClick},
  {id: 'u2demo-inspector/2 an entity chip click shows the entity, not the chip control', run: checkChipShowsEntity},
  {id: 'u2demo-inspector/3 clicking bare whitespace leaves the panel as it is', run: checkWhitespace},
  {id: 'u2demo-inspector/4 an unregistered control says it has no properties', run: checkUnregisteredControl},
];
