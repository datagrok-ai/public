/* Locators and gestures against a static page: the kinds and the platform names as the
   bindings register them, resolved on markup shaped like the u2 and Dart contracts. Runs in the
   library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

/* These tests need a browser, and the other test files do not: the launch is guarded, so a
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
import {noSpaceOnServer} from '../bindings/platform/steps.js';
import {el} from '../src/runtime/args.js';
import {clear, keysOf, press, pressIn, readExpanded, select, setExpanded, typeInto, withKeys} from '../src/runtime/gestures.js';
import {atFeatureEnd, feature, takeErrors, watchErrors} from '../src/runtime/harness.js';
import {expectState, expectText} from '../src/runtime/assertions.js';
import {whileExpectedToFail} from '../src/runtime/patience.js';
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
<div class="d4-toolbox"><div class="grok-toolbox-header">Pin toolbox</div></div>
<div class="d4-toolbox" caption=" ">
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
<div class="panel-base">
  <div class="panel-titlebar">
    <div class="panel-titlebar-text">Filters</div>
    <i class="grok-icon grok-font-icon-help" name="icon-font-icon-help">titlebar</i>
    <i class="grok-icon grok-font-icon-settings" name="icon-font-icon-settings"></i>
  </div>
  <div name="viewer-Filters">
    <div class="d4-filter-group-header"><i class="grok-icon fal fa-question" name="icon-question">group</i></div>
  </div>
</div>
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

const PANE = (button: boolean) => '<div class="d4-accordion-pane" name="pane-Grants">' +
  '<div class="d4-accordion-pane-header" name="div-section--Grants">Grants</div>' +
  `<div class="d4-accordion-pane-content">${button ? '<button class="ui-btn">MANAGE</button>' : ''}</div></div>`;

async function withHost(html: string, body: () => Promise<void>): Promise<void> {
  await page!.evaluate((h) => {
    const host = document.createElement('div');
    host.id = 'late-host';
    host.innerHTML = h;
    document.body.appendChild(host);
  }, html);
  try {
    await body();
  }
  finally {
    await page!.evaluate(() => document.getElementById('late-host')?.remove());
  }
}

const setHost = (html: string) => page!.evaluate((h) => { document.getElementById('late-host')!.innerHTML = h; }, html);

scenario('a scope read while its container rebuilds finds the target once the container is back', () =>
  withHost('<div class="d4-accordion-pane-header" name="div-section--Grants">Grants</div>', async () => {
    const button = await locateActionable(page!, el('MANAGE button in "Grants" section'));
    await setHost(PANE(true));
    assert.equal(await button.count(), 1);
  }));

scenario('a target that renders late in a whole scope, and in an ordinal one, is found in that scope', () =>
  withHost(PANE(false) + PANE(false), async () => {
    const inFirst = await locateActionable(page!, el('MANAGE button in "Grants" section'));
    const inSecond = await locateActionable(page!, el('MANAGE button in second "Grants" section'));
    await setHost(PANE(false) + PANE(true));
    assert.equal(await inSecond.count(), 1);
    assert.equal(await inFirst.count(), 1);
    assert.equal(await page!.locator('#late-host [name="pane-Grants"]').nth(1).locator('button').count(), 1);
  }));

scenario('ordinals count the visible matches, as a gesture does', async () => {
  assert.equal(await text('second item in results list'), 'beta');
  assert.equal(await text('last item in results list'), 'beta');
  assert.equal(await (await locateActionable(page!, el('item in results list'))).count(), 2);
});

scenario('a menu item matches its own label, not its children', async () => {
  assert.equal(await count('"As CSV" menu item in context menu'), 1);
  assert.equal(await (await locate(page!, el('Export menu item in context menu'))).getAttribute('name'), 'div-Export');
});

scenario('the help icon is the one of the title bar around the viewer, not a "?" inside it', async () => {
  assert.equal(await count('help icon of filter panel'), 1);
  assert.equal(await text('help icon of filter panel'), 'titlebar');
  assert.equal(await text('help icon of Filters viewer'), 'titlebar');
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

scenario('typing replaces existing text in a focused editor on every platform', async () => {
  const input = page!.locator('[data-u2-name="org"] input');
  await input.fill('45');
  await input.press('End');
  await typeInto(page!, el('org input'), '99');
  assert.equal(await input.inputValue(), '99');
  await typeInto(page!, el('org input'), '7');
  assert.equal(await input.inputValue(), '7');
});

scenario('clearing removes all existing text on every platform', async () => {
  const input = page!.locator('[data-u2-name="org"] input');
  await input.fill('0.00');
  await clear(page!, el('org input'));
  assert.equal(await input.inputValue(), '');
});

scenario('the portable selection modifier clicks without a context menu and is released', async () => {
  const button = page!.getByRole('button', {name: 'OK', exact: true});
  await button.evaluate((el) => {
    for (const type of ['click', 'contextmenu'])
      el.addEventListener(type, (event) => {
        const e = event as MouseEvent;
        el.setAttribute('data-gesture', `${e.type}:${e.ctrlKey || e.metaKey}:${e.shiftKey}`);
        e.preventDefault();
      });
  });
  for (const keys of [keysOf('ctrl+shift'), ['Control', 'Shift'], ['ControlOrMeta', 'Shift']]) {
    await withKeys(page!, keys, () => button.click());
    assert.equal(await button.getAttribute('data-gesture'), 'click:true:true');
  }
  await button.click();
  assert.equal(await button.getAttribute('data-gesture'), 'click:false:false');
});

scenario('Control shortcuts select text while ControlLeft sends physical Control', async () => {
  const input = page!.locator('[data-u2-name="org"] input');
  await input.fill('abcdef');
  await input.press('End');
  await press(page!, 'Control+A');
  assert.equal(await input.evaluate((el: HTMLInputElement) => el.selectionEnd! - el.selectionStart!), 6);
  await input.evaluate((el) => el.addEventListener('keydown', (event) => {
    const e = event as KeyboardEvent;
    el.setAttribute('data-modifiers', `${e.ctrlKey}:${e.metaKey}:${e.shiftKey}`);
  }));
  await pressIn(page!, el('org input'), 'ControlLeft+Shift+C');
  assert.equal(await input.getAttribute('data-modifiers'), 'true:false:true');
});

scenario('Delete and Del follow Datagrok commands while ForwardDelete and Backspace stay literal', async () => {
  const input = page!.locator('[data-u2-name="org"] input');
  await input.evaluate((el) => {
    el.addEventListener('keydown', (event) => {
      const e = event as KeyboardEvent;
      el.setAttribute('data-key', `${e.key}:${e.keyCode}:${e.shiftKey}`);
    });
  });
  const mac = await page!.evaluate(() => navigator.platform.startsWith('Mac'));
  await pressIn(page!, el('org input'), 'Shift+Del');
  assert.equal(await input.getAttribute('data-key'), mac ? 'Backspace:8:true' : 'Delete:46:true');
  await pressIn(page!, el('org input'), 'Shift+Delete');
  assert.equal(await input.getAttribute('data-key'), mac ? 'Backspace:8:true' : 'Delete:46:true');
  await pressIn(page!, el('org input'), 'Shift+ForwardDelete');
  assert.equal(await input.getAttribute('data-key'), 'Delete:46:true');
  await pressIn(page!, el('org input'), 'Shift+Backspace');
  assert.equal(await input.getAttribute('data-key'), 'Backspace:8:true');
});


scenario('select handles both a named native select and a host with a select inside', async () => {
  await page!.setContent(`<select class="ui-input-root" name="input-host-Engine"><option>A</option><option>B</option></select>
    <div class="ui-input-root" name="input-host-Mode"><select><option>C</option><option>D</option></select></div>`);
  await select(page!, el('Engine input'), 'B');
  await select(page!, el('Mode input'), 'D');
  assert.deepEqual(await page!.locator('select').evaluateAll((all) => all.map((e) => (e as HTMLSelectElement).value)), ['B', 'D']);
});

scenario('property category expansion and collapse are idempotent', async () => {
  await page!.setContent(`<table><tr class="property-grid-item property-grid-category" name="prop-category-Axes"><td><i class="property-grid-icon-plus"></i><span class="property-grid-item-name-text">Axes</span></td></tr></table>`);
  await page!.locator('.property-grid-category').evaluate((e) => e.addEventListener('click', () => {
    const icon = e.querySelector('i')!;
    icon.classList.toggle('property-grid-icon-minus');
    icon.classList.toggle('property-grid-icon-plus');
    e.setAttribute('data-clicks', String(Number(e.getAttribute('data-clicks')) + 1));
  }));
  for (const expanded of [true, true, false, false]) {
    await setExpanded(page!, el('Axes category'), expanded);
    assert.equal(await readExpanded(page!.locator('.property-grid-category')), expanded);
  }
  assert.equal(await page!.locator('.property-grid-category').getAttribute('data-clicks'), '2');
});

scenario('ready requires an explicit completed and valid asynchronous result', async () => {
  await page!.setContent('<div class="d4-pm-view-preview">preview</div>');
  const preview = page!.locator('.d4-pm-view-preview');
  await expectState(page!, el('model preview'), 'ready', true);
  await preview.evaluate((e) => e.setAttribute('aria-busy', 'true'));
  await expectState(page!, el('model preview'), 'ready', true);
  await preview.evaluate((e) => { e.setAttribute('aria-busy', 'false'); e.setAttribute('aria-invalid', 'true'); });
  await expectState(page!, el('model preview'), 'ready', true);
  await preview.evaluate((e) => e.setAttribute('aria-invalid', 'false'));
  await expectState(page!, el('model preview'), 'ready');
});

scenario('tooltip text assertions ignore hidden content retained between hovers', async () => {
  await page!.setContent('<div class="d4-tooltip" style="display:none">Pearson R: 0.4</div>');
  await expectText(page!, el('tooltip'), 'Pearson R', {negate: true});
  await assert.rejects(whileExpectedToFail(() => expectText(page!, el('tooltip'), 'Pearson R')));
  await page!.locator('.d4-tooltip').evaluate((e) => (e as HTMLElement).style.display = 'block');
  await expectText(page!, el('tooltip'), 'Pearson R: 0.4', {exact: true});
  await assert.rejects(whileExpectedToFail(() => expectText(page!, el('tooltip'), 'Pearson R', {negate: true})));
  await page!.locator('.d4-tooltip').evaluate((e) => e.remove());
  await expectText(page!, el('tooltip'), 'Pearson R', {negate: true});
});

scenario('space cleanup deletes existing fixtures even when server filters return no matches', async () => {
  let afterAll!: () => Promise<void>;
  const api = {afterEach: () => undefined, afterAll: (hook: () => Promise<void>) => { afterAll = hook; }};
  const session = feature(api as unknown as Parameters<typeof feature>[0]);
  const cleanupPage = await session.page(browser!);
  await cleanupPage.setContent('<div class="grok-browse-icons"><i class="fa fa-sync" title="Refresh">Refresh</i></div>' +
    '<div class="layout-browse"><span id="stale-fixture">BDD Fixture</span></div>');
  await cleanupPage.locator('[title="Refresh"]').evaluate((icon) => {
    icon.addEventListener('click', () => document.getElementById('stale-fixture')!.remove());
  });
  await cleanupPage.evaluate(() => {
    let spaces = [
      {id: 'fixture', name: 'BDDFixture', friendlyName: 'BDD Fixture'},
      {id: 'unrelated', name: 'Keep', friendlyName: 'Keep'},
    ];
    const data = {
      order() { return data; },
      filter() { return {async list() { return []; }}; },
      async list() { return spaces; },
      async find(id: string) { return spaces.find((space) => space.id === id); },
      async delete(space: {id: string}) { spaces = spaces.filter((item) => item.id !== space.id); },
      async createRootSpace(name: string) { spaces.push({id: 'new-fixture', name, friendlyName: name}); },
    };
    (window as any).grok = {dapi: {spaces: data}};
  });
  await noSpaceOnServer(cleanupPage, 'BDD Fixture');
  const remaining = () => cleanupPage.evaluate(async () =>
    (await (window as any).grok.dapi.spaces.list()).map((space: {id: string}) => space.id));
  assert.deepEqual(await remaining(), ['unrelated']);
  assert.equal(await cleanupPage.locator('#stale-fixture').count(), 0);
  await cleanupPage.evaluate(() => (window as any).grok.dapi.spaces.createRootSpace('BDD Fixture'));
  assert.deepEqual(await remaining(), ['unrelated', 'new-fixture']);
  await afterAll();
  assert.deepEqual(await remaining(), ['unrelated']);
});

scenario('feature teardown attempts every cleanup and reports synchronous and asynchronous failures', async () => {
  let afterAll!: () => Promise<void>;
  const api = {afterEach: () => undefined, afterAll: (hook: () => Promise<void>) => { afterAll = hook; }};
  const session = feature(api as unknown as Parameters<typeof feature>[0]);
  const cleanupPage = await session.page(browser!);
  const ran: number[] = [];
  const errors = [new Error('synchronous'), new Error('asynchronous')];
  atFeatureEnd(cleanupPage, () => { ran.push(1); throw errors[0]; });
  atFeatureEnd(cleanupPage, async () => { ran.push(2); throw errors[1]; });
  atFeatureEnd(cleanupPage, async () => { ran.push(3); });
  await assert.rejects(afterAll, (e: unknown) => e instanceof AggregateError &&
    e.errors[0] === errors[0] && e.errors[1] === errors[1] && e.errors.length === 2);
  assert.deepEqual(ran, [1, 2, 3]);
  await afterAll();
});

test('a run suffix is stable within a feature and distinct across feature instances', () => {
  const api = {afterEach: () => undefined, afterAll: () => undefined};
  const a = feature(api as unknown as Parameters<typeof feature>[0]);
  const b = feature(api as unknown as Parameters<typeof feature>[0]);
  assert.equal(a.text('literal'), 'literal');
  assert.equal(a.text('{run}/{run}').split('/')[0], a.text('{run}'));
  assert.equal(a.text('{run}/{run}').split('/')[1], a.text('{run}'));
  assert.notEqual(a.text('{run}'), b.text('{run}'));
  assert.match(a.text('{run}'), /^[0-9a-f-]{36}$/);
  assert.equal(a.text('u{time}.{time}').split('.')[0], `u${a.text('{time}')}`);
  assert.match(a.text('{time}'), /^\d{13,}$/);
  assert.notEqual(a.text('{time}'), b.text('{time}'));
});

test('the translated stack trace of a logged error is not a second error', () => {
  const listeners: ((m: unknown) => void)[] = [];
  const fake = {on: (event: string, fn: (m: unknown) => void) => { if (event === 'console') listeners.push(fn); }} as unknown as Page;
  const log = (text: string) => listeners.forEach((fn) => fn({type: () => 'error', text: () => text, location: () => ({url: ''})}));
  watchErrors(fake);
  log('NullError: y on null\nTranslating stack trace... Look below, ID = Cd2');
  log('Stack trace Cd2\n\tpackages/d4/src/w.dart 1:1');
  const joined = takeErrors(fake);
  assert.equal(joined.length, 1);
  assert.match(joined[0], /w\.dart/);
  log('NullError: x on null\nTranslating stack trace... Look below, ID = Ab1');
  assert.equal(takeErrors(fake).length, 1);
  log('Stack trace Ab1\n\tpackages/d4/src/x.dart 1:1');
  assert.deepEqual(takeErrors(fake), []);
  log('Stack trace Zz9\n\tpackages/d4/src/y.dart 1:1');
  assert.equal(takeErrors(fake).length, 1);
});
