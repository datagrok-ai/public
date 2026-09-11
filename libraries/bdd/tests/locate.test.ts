/* The phrase-to-locator layer against a static page: the kinds and the platform names as the
   bindings register them, resolved on markup shaped like the u2 and Dart contracts. Runs in the
   library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

/* These eight need a browser, and the other five test files do not: the launch is guarded, so a
   runner without Playwright's Chromium (the libraries CI installs none) skips them and says why,
   and the rest of the unit tests still run. */
let missing = '';
const scenario = (name: string, fn: () => Promise<void>): void => void test(name, async (t) => {
  if (page === undefined) {
    t.skip(missing || 'no page');
    return;
  }
  await fn();
});
import '../bindings/common/kinds.js';
import '../bindings/platform/elements.js';
import {el} from '../src/runtime/args.js';
import {locate, locateActionable} from '../src/runtime/locate.js';

const PAGE = `
<div data-u2="dialog" data-u2-name="loadProjectDialog">
  <div class="u2-dialog-title"><span class="u2-dialog-title-text">Load Project</span></div>
  <div data-u2="form" data-u2-name="project">
    <div data-u2="text-input" data-u2-name="org">
      <label data-u2-part="label">Organization</label>
      <span data-u2-part="editor"><input value="acme"></span>
    </div>
    <div data-u2="text-input">
      <label data-u2-part="label">First name</label>
      <span data-u2-part="editor"><input value="Ada"></span>
    </div>
    <div data-u2="choice-input" data-u2-name="scope"><select><option>Mine</option><option selected>All</option></select></div>
  </div>
  <button>OK</button><button>CANCEL</button>
</div>
<div class="d4-toolbox">
  <div name="div-section--Viewers">
    <div name="icon-scatter-plot"></div><div name="icon-bar-chart"></div>
  </div>
</div>
<div class="ui-input-root ui-input-text" name="input-host-Caption"><label class="ui-input-label">Caption</label><input class="ui-input-editor" value="x"></div>
<div data-u2="list" data-u2-name="results">
  <div class="u2-list-row" role="option">alpha</div>
  <div class="u2-list-row" role="option" aria-selected="true">beta</div>
  <div class="u2-list-row" role="option" style="display:none">gamma</div>
</div>
<div class="d4-menu-popup">
  <div class="d4-menu-item" name="div-Export" role="menuitem"><div class="d4-menu-item-label">Export</div>
    <div class="d4-menu-item" name="div-Export---As-CSV" role="menuitem"><div class="d4-menu-item-label">As CSV</div></div>
  </div>
</div>
<div data-u2="menu" data-u2-owner="scope"><div class="u2-menu-item" role="menuitem"><span class="u2-menu-label">Everything</span></div></div>
`;

let browser: Browser | undefined;
let page: Page | undefined;

before(async () => {
  try {
    browser = await chromium.launch();
  }
  catch (e) {
    missing = `Chromium does not start here: ${String(e).slice(0, 80)}`;
    return;
  }
  page = await browser.newPage();
  await page.setContent(PAGE);
});

after(async () => {
  await browser?.close();
});

const count = async (phrase: string): Promise<number> => (await locate(page!, el(phrase))).count();
const text = async (phrase: string): Promise<string> => (await (await locate(page!, el(phrase))).first().textContent()) ?? '';

scenario('a kind by its data-u2-name, and by its label', async () => {
  assert.equal(await count('org input'), 1);
  assert.equal(await count('organization text input'), 1);
  assert.equal(await count('"First name" input'), 1);
});

scenario('a part of an input, and a part of a dialog', async () => {
  assert.equal(await (await locate(page!, el('editor of org input'))).locator('input').inputValue(), 'acme');
  assert.equal(await text('title of load project dialog'), 'Load Project');
});

scenario('composition scopes the inner phrase inside the outer one', async () => {
  assert.equal(await count('org input in project form in load project dialog'), 1);
  assert.equal(await count('OK button in load project dialog'), 1);
  assert.equal(await count('OK button in project form'), 0);
});

scenario('the Dart name conventions and the platform names', async () => {
  assert.equal(await count('scatter plot icon on toolbox'), 1);
  assert.equal(await count('scatter plot icon in viewers section of toolbox'), 1);
  assert.equal(await count('Caption input'), 1);
});

scenario('ordinals, and the visible matches a gesture acts on', async () => {
  assert.equal(await text('second item in results list'), 'beta');
  assert.equal(await text('last item in results list'), 'gamma');
  assert.equal(await (await locateActionable(page!, el('item in results list'))).count(), 2);
});

scenario('a menu item matches its own label, not its children', async () => {
  assert.equal(await count('"As CSV" menu item in context menu'), 1);
  assert.equal(await (await locate(page!, el('Export menu item in context menu'))).getAttribute('name'), 'div-Export');
});

scenario('a popup portaled out of its owner is found through the owner edge', async () => {
  assert.equal(await count('Everything menu item in scope input'), 1);
});

scenario('a phrase that matches nothing still has a locator to fail against, at once', async () => {
  const start = Date.now();
  assert.equal(await count('nowhere button in load project dialog'), 0);
  assert.equal(await count('OK button in missing dialog'), 0);
  assert.ok(Date.now() - start < 5000, 'a scope that is not there must not be waited for');
});
