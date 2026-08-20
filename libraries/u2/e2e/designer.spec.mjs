/* The u2 designer live-pass checklist, as deterministic checks against local mode. Nodes are
   addressed by `data-u2-name` (WO-T1) wherever the spec names them; the platform chrome is addressed
   by its own classes. Run through run.mjs — `npm run e2e:local`. */
import {ok, shot} from './local.mjs';

/** What U2Demo's `u2DesignerApp` opens on (packages/U2Demo/src/designer.ts). */
export const APP = {package: 'U2demo', func: 'u2DesignerApp'};

const ROW = '.u2-designer-toolbox .u2-tree-row';
const CAPTION = '.u2-designer-caption';
/** The status bar panel carries no class of its own — it is the leaf div reporting the node count. */
const statusText = (page) => page.evaluate(() => [...document.querySelectorAll('div')]
  .find((x) => x.children.length === 0 && /\d+ node/.test(x.textContent))?.textContent ?? '');

const caption = (page) => page.evaluate((sel) =>
  document.querySelector(sel)?.textContent ?? '', CAPTION);

const balloons = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .map((b) => b.textContent.trim()).filter((t) => t));

const clearBalloons = (page) => page.evaluate(() =>
  document.querySelectorAll('.d4-balloon-content, .d4-balloon').forEach((b) => b.remove()));

/** The property panel as `{section: {row: value}}` — what the context panel actually shows. */
const panel = (page) => page.evaluate(() => {
  const props = document.querySelector('.u2-designer-properties');
  const sections = {};
  let current = null;
  for (const el of props?.children ?? []) {
    if (el.tagName === 'H3')
      sections[current = el.textContent] = {};
    else if (current) {
      for (const tr of el.querySelectorAll('tr'))
        sections[current][tr.children[0]?.textContent.trim()] = tr.children[1]?.textContent.trim();
    }
  }
  return sections;
});

/** The panel is asynchronous — the shell renders it after `grok.shell.o` is set. */
async function waitCaption(page, name, timeout = 8000) {
  await page.waitForFunction(({sel, name}) =>
    (document.querySelector(sel)?.textContent ?? '').startsWith(`${name} (`),
  {sel: CAPTION, name}, {timeout}).catch(() => {});
  return caption(page);
}

/** The designer's own status bar: written synchronously by the selection, so it is the authority
 * on what the designer selected — the panel is the platform rendering that, one step behind. */
async function waitStatus(page, name, timeout = 3000) {
  await page.waitForFunction((name) => [...document.querySelectorAll('div')]
    .some((x) => x.children.length === 0 && new RegExp(`${name} · \\d+ node`).test(x.textContent)),
  name, {timeout}).catch(() => {});
  return statusText(page);
}

/** One twistie at a time: expanding re-renders the virtual list, which invalidates every other
 * handle held across the click. */
async function expandTree(page) {
  const collapsed = page.locator('.u2-designer-toolbox .u2-tree-twistie:not(.u2-tree-twistie-hidden)')
    .filter({hasText: '▸'});
  for (let i = 0; i < 30 && await collapsed.count() > 0; i++) {
    await collapsed.first().click({timeout: 2000}).catch(() => {});
    await page.waitForTimeout(80);
  }
}

const row = (page, name) => page.locator(ROW).filter({has: page.locator('.u2-tree-label', {hasText: name})}).first();

/** Selects a node from the structure tree. The platform's context panel runs one `grok.shell.o`
 * assignment behind (verified: reading `grok.shell.o` right after a single write still answers the
 * previous object), so a selection that has reached the status bar but not the panel is re-asserted
 * with the click a user would make — the recovery check below covers that path on its own. */
async function selectRow(page, name) {
  await row(page, name).click();
  await waitStatus(page, name);
  const shown = await waitCaption(page, name, 1500);
  if (shown.startsWith(`${name} (`))
    return shown;
  await row(page, name).click();
  return waitCaption(page, name);
}

async function checkView(page) {
  await expandTree(page);
  await shot(page, '1-view');
  const view = await page.evaluate(() => {
    const surface = document.querySelector('.u2-designer-surface');
    const tree = document.querySelector('.u2-designer-toolbox .u2-tree');
    return {
      designer: !!document.querySelector('.u2-designer'),
      surfaceKids: surface ? surface.children.length : -1,
      treeHeight: tree ? Math.round(tree.getBoundingClientRect().height) : -1,
      labels: [...document.querySelectorAll('.u2-designer-toolbox .u2-tree-label')].map((l) => l.textContent),
      named: [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name),
      ribbon: [...document.querySelectorAll('.d4-ribbon button, .d4-ribbon .u2-btn')]
        .map((b) => b.textContent.trim()).filter((t) => t),
    };
  });
  const status = await statusText(page);
  ok('1a/canvas-renders-spec', view.designer && view.surfaceKids > 0, `surface children=${view.surfaceKids}`);
  ok('1b/structure-tree', view.treeHeight > 80 && view.labels.includes('saveButton'),
    `tree=${view.treeHeight}px labels=${JSON.stringify(view.labels)}`);
  ok('1c/ribbon-actions', ['Open…', 'Copy spec', 'Design', 'Run'].every((t) => view.ribbon.includes(t)),
    view.ribbon.join('|'));
  ok('1d/status-path-and-count', /^layout · \d+ nodes/.test(status), status);
  ok('1e/named-nodes-stamped', ['nameInput', 'doseInput', 'saveButton'].every((n) => view.named.includes(n)),
    JSON.stringify(view.named));
}

async function checkPanel(page) {
  const buttonCaption = await selectRow(page, 'saveButton');
  const button = await panel(page);
  ok('2a/panel-identity', buttonCaption === 'saveButton (u2-button)', `caption="${buttonCaption}"`);
  ok('2b/props-carry-real-values', button.Properties?.text === 'Save',
    JSON.stringify(button.Properties ?? null));
  ok('2c/events-listed', button.Events?.click === 'cmd:save', JSON.stringify(button.Events ?? null));

  await selectRow(page, 'layout');
  const splitter = await panel(page);
  ok('2d/props-read-through-spec-props', splitter.Properties?.direction === 'horizontal',
    JSON.stringify(splitter.Properties ?? null));

  const boundCaption = await selectRow(page, 'nameInput');
  const bound = await panel(page);
  await shot(page, '2-panel');
  ok('2e/bindings-section', bound.Bindings?.value === '$.reagent',
    `caption="${boundCaption}" bindings=${JSON.stringify(bound.Bindings ?? null)}`);
}

/** The click goes through the glass pane, the form and the input wrapper to the innermost node the
 * spec owns — `elementsFromPoint` filtered by the canvas, not each control's own handlers. */
async function checkHitTest(page) {
  const box = await page.locator('.u2-designer-surface [data-u2-name="doseInput"] input').first().boundingBox();
  if (!box) {
    ok('3a/innermost-node-selected', false, 'doseInput is not on the canvas');
    return;
  }
  const click = () => page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await click();
  const status = await waitStatus(page, 'doseInput');
  let hit = await waitCaption(page, 'doseInput', 1500);
  if (!hit.startsWith('doseInput (')) {
    await click();
    hit = await waitCaption(page, 'doseInput');
  }
  await shot(page, '3-hit-test');
  ok('3a/innermost-node-selected', hit === 'doseInput (u2-number-input)', `caption="${hit}"`);
  ok('3b/status-follows-hit', /› doseInput · \d+ nodes/.test(status), status);
  ok('3c/adorner-visible', await page.evaluate(() =>
    getComputedStyle(document.querySelector('.u2-designer-selected')).display) === 'block');
}

async function checkAdorner(page) {
  const geometry = () => page.evaluate(() => {
    const round = (r) => r && {x: Math.round(r.x), y: Math.round(r.y),
      w: Math.round(r.width), h: Math.round(r.height)};
    return {
      box: round(document.querySelector('.u2-designer-selected').getBoundingClientRect()),
      target: round(document.querySelector('.u2-designer-surface [data-u2-name="doseInput"]')
        ?.getBoundingClientRect()),
    };
  });
  const near = (g) => !!g.target && ['x', 'y', 'w', 'h']
    .every((k) => Math.abs(g.box[k] - g.target[k]) <= 3);
  const before = await geometry();
  ok('4a/adorner-matches-selection', near(before), JSON.stringify(before));

  await page.setViewportSize({width: 1200, height: 780});
  await page.waitForTimeout(400);
  const resized = await geometry();
  ok('4b/adorner-follows-resize', near(resized), JSON.stringify(resized));
  await shot(page, '4-adorner');
  await page.setViewportSize({width: 1600, height: 1000});
  await page.waitForTimeout(400);
}

async function checkModes(page) {
  const input = page.locator('.u2-designer-surface [data-u2-name="nameInput"] input').first();
  const before = await input.inputValue();
  const box = await input.boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.keyboard.type('ZZZ', {delay: 20});
  ok('5a/design-mode-inert', await input.inputValue() === before,
    `before="${before}" after="${await input.inputValue()}"`);

  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(300);
  const running = await page.evaluate(() => ({
    pane: !!document.querySelector('.u2-designer-pane'),
    running: document.querySelector('.u2-designer').classList.contains('u2-designer-running'),
    adorner: getComputedStyle(document.querySelector('.u2-designer-selected')).display,
  }));
  ok('5b/run-mode-drops-glass-pane', !running.pane && running.running && running.adorner === 'none',
    JSON.stringify(running));

  await input.click();
  await input.press('End');
  await page.keyboard.type('ZZZ', {delay: 20});
  const typed = await input.inputValue();
  ok('5c/run-mode-accepts-typing', typed !== before && /ZZZ/.test(typed), `"${before}" -> "${typed}"`);

  await clearBalloons(page);
  await page.locator('.u2-designer-surface [data-u2-name="saveButton"]').first().click();
  await page.waitForFunction(() => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => /cmd:save ran/.test(b.textContent)), null, {timeout: 5000}).catch(() => {});
  const fired = await balloons(page);
  await shot(page, '5-run-mode');
  ok('5d/command-fires-in-run-mode', fired.some((b) => /cmd:save ran/.test(b)), JSON.stringify(fired));

  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
  ok('5e/back-to-design', await page.evaluate(() => !!document.querySelector('.u2-designer-pane')));
}

/** The 2b fix: a click on the row that is already selected moves no signal, so the view has to
 * re-assert the selection itself — otherwise the panel keeps whatever it showed before. */
async function checkReclick(page) {
  await selectRow(page, 'nameInput');
  await selectRow(page, 'notesInput');
  await page.evaluate(() => {
    grok.shell.o = grok.shell.user;
  });
  await page.waitForTimeout(300);
  await row(page, 'notesInput').click();
  const back = await waitCaption(page, 'notesInput');
  await shot(page, '6-reclick');
  ok('6a/same-row-reclick-recovers', back === 'notesInput (u2-text-area)', `caption="${back}"`);
  ok('6b/status-agrees', /notesInput/.test(await statusText(page)), await statusText(page));
}

const ROUND_TRIP = JSON.stringify({
  $schema: 'dg-ui/1',
  root: {tag: 'u2-panel', name: 'roundTrip', children: [
    {tag: 'h2', props: {text: 'Round trip'}},
    {tag: 'u2-form', name: 'rtForm', children: [
      {tag: 'u2-text-input', name: 'rtName', props: {label: 'Reagent', value: 'Ethanol'}},
      {tag: 'u2-bool-input', name: 'rtFlag', props: {label: 'Checked', value: false}},
    ]},
  ]},
}, null, 2);

async function checkOpen(page) {
  await page.locator('.d4-ribbon').getByText('Open', {exact: false}).first().click();
  await page.waitForSelector('.u2-designer-spec-editor textarea', {timeout: 10000});
  const prefilled = await page.locator('.u2-designer-spec-editor textarea').inputValue();
  ok('7a/open-dialog-prefilled', prefilled.includes('"$schema": "dg-ui/1"') && prefilled.includes('layout'),
    `${prefilled.length} chars`);

  await page.locator('.u2-designer-spec-editor textarea').fill(ROUND_TRIP);
  await page.locator('.d4-dialog button, .u2-dialog button').filter({hasText: /^OK$/i}).first().click();
  // the selection is read before the tree is expanded: expanding a row selects it
  const status = await waitStatus(page, 'roundTrip');
  await expandTree(page);
  await shot(page, '7-round-trip');
  const after = await page.evaluate(() => ({
    named: [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name),
    value: document.querySelector('.u2-designer-surface [data-u2-name="rtName"] input')?.value,
    labels: [...document.querySelectorAll('.u2-designer-toolbox .u2-tree-label')].map((l) => l.textContent),
  }));
  ok('7b/canvas-rerenders', after.value === 'Ethanol' && after.named.includes('rtFlag'),
    `value=${after.value} named=${JSON.stringify(after.named)}`);
  ok('7c/tree-and-status-follow', after.labels.includes('rtForm') && /^roundTrip.* · 5 nodes/.test(status),
    `labels=${JSON.stringify(after.labels)} status="${status}"`);
  // the status bar is u2's authority on the selection (the platform panel renders one
  // current-object change behind and a lone programmatic set never renders — a known
  // platform defect; the click below is the user gesture that converges it, as elsewhere)
  const opened = await selectRow(page, 'roundTrip');
  ok('7d/root-selected', /^roundTrip · /.test(status) && opened === 'roundTrip (u2-panel)',
    `status="${status}" caption-after-click="${opened}"`);

  await clearBalloons(page);
  await page.locator('.d4-ribbon').getByText('Open', {exact: false}).first().click();
  await page.waitForSelector('.u2-designer-spec-editor textarea', {timeout: 10000});
  await page.locator('.u2-designer-spec-editor textarea').fill('{"nope": 1}');
  await page.locator('.d4-dialog button, .u2-dialog button').filter({hasText: /^OK$/i}).first().click();
  await page.waitForTimeout(600);
  const refused = await balloons(page);
  const kept = await page.evaluate(() =>
    !!document.querySelector('.u2-designer-surface [data-u2-name="rtName"]'));
  ok('7e/invalid-spec-refused', kept && refused.some((b) => /\$schema/.test(b)),
    `canvas kept=${kept} balloons=${JSON.stringify(refused)}`);
}

export const checks = [
  {id: '1 view renders', run: checkView},
  {id: '2 selection to context panel', run: checkPanel},
  {id: '3 canvas hit-test', run: checkHitTest},
  {id: '4 adorner geometry', run: checkAdorner},
  {id: '5 design inertness and run liveness', run: checkModes},
  {id: '6 same-row re-click', run: checkReclick},
  {id: '7 open round-trip', run: checkOpen},
];
