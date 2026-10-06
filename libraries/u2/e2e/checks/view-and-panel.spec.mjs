/* The live-pass checklist on the spec the host wrote: the canvas, the structure tree, the context
   panel, the hit-test, the adorner, the two modes, the same-row re-click and the Open round trip. */
import {ok, shot} from '../local.mjs';
import {ROUND_TRIP, balloons, clearBalloons, expandTree, openSpec, panel, reopenApp, row, selectRow,
  statusText, waitCaption, waitStatus} from '../lib.mjs';

/** The spec `u2DesignerApp` opens on, as the host hands it over — not whatever the previous file left. */
export async function fixture(page) {
  await reopenApp(page);
}

async function checkView(page) {
  await expandTree(page);
  await shot(page, 'view-and-panel-1-view');
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
      // the icon ribbon (F15): the verbs are addressed by accessible name, not text
      ribbonNames: [...document.querySelectorAll('.d4-ribbon button')]
        .map((b) => b.getAttribute('aria-label') ?? '').filter((t) => t),
    };
  });
  const status = await statusText(page);
  ok('view-and-panel/1a/canvas-renders-spec', view.designer && view.surfaceKids > 0,
    `surface children=${view.surfaceKids}`);
  ok('view-and-panel/1b/structure-tree', view.treeHeight > 80 && view.labels.includes('saveButton'),
    `tree=${view.treeHeight}px labels=${JSON.stringify(view.labels)}`);
  ok('view-and-panel/1c/ribbon-actions', ['Open', 'Copy spec'].every((t) => view.ribbonNames.includes(t)) &&
    ['Design', 'Run'].every((t) => view.ribbon.includes(t)),
  `names=${view.ribbonNames.join('|')} text=${view.ribbon.join('|')}`);
  ok('view-and-panel/1d/status-path-and-count', /^layout · \d+ nodes/.test(status), status);
  ok('view-and-panel/1e/named-nodes-stamped',
    ['nameInput', 'doseInput', 'saveButton'].every((n) => view.named.includes(n)),
    JSON.stringify(view.named));
}

async function checkPanel(page) {
  const buttonCaption = await selectRow(page, 'saveButton');
  const button = await panel(page);
  ok('view-and-panel/2a/panel-identity', buttonCaption === 'saveButton (u2-button)', `caption="${buttonCaption}"`);
  ok('view-and-panel/2b/props-carry-real-values', button.Properties?.text === 'Save',
    JSON.stringify(button.Properties ?? null));
  ok('view-and-panel/2c/events-listed', button.Events?.click === 'cmd:save', JSON.stringify(button.Events ?? null));

  await selectRow(page, 'layout');
  const splitter = await panel(page);
  ok('view-and-panel/2d/props-read-through-spec-props', splitter.Properties?.direction === 'horizontal',
    JSON.stringify(splitter.Properties ?? null));

  const boundCaption = await selectRow(page, 'nameInput');
  const bound = await panel(page);
  await shot(page, 'view-and-panel-2-panel');
  ok('view-and-panel/2e/bindings-section', bound.Bindings?.value === '$.reagent',
    `caption="${boundCaption}" bindings=${JSON.stringify(bound.Bindings ?? null)}`);
}

/** The click goes through the glass pane, the form and the input wrapper to the innermost node the
 * spec owns — `elementsFromPoint` filtered by the canvas, not each control's own handlers. */
async function checkHitTest(page) {
  const box = await page.locator('.u2-designer-surface [data-u2-name="doseInput"] input').first().boundingBox();
  if (!box) {
    ok('view-and-panel/3a/innermost-node-selected', false, 'doseInput is not on the canvas');
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
  await shot(page, 'view-and-panel-3-hit-test');
  ok('view-and-panel/3a/innermost-node-selected', hit === 'doseInput (u2-number-input)', `caption="${hit}"`);
  ok('view-and-panel/3b/status-follows-hit', /› doseInput · \d+ nodes/.test(status), status);
  ok('view-and-panel/3c/adorner-visible', await page.evaluate(() =>
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
  ok('view-and-panel/4a/adorner-matches-selection', near(before), JSON.stringify(before));

  await page.setViewportSize({width: 1200, height: 780});
  await page.waitForTimeout(400);
  const resized = await geometry();
  ok('view-and-panel/4b/adorner-follows-resize', near(resized), JSON.stringify(resized));
  await shot(page, 'view-and-panel-4-adorner');
  await page.setViewportSize({width: 1600, height: 1000});
  await page.waitForTimeout(400);
}

async function checkModes(page) {
  const input = page.locator('.u2-designer-surface [data-u2-name="nameInput"] input').first();
  const before = await input.inputValue();
  const box = await input.boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.keyboard.type('ZZZ', {delay: 20});
  ok('view-and-panel/5a/design-mode-inert', await input.inputValue() === before,
    `before="${before}" after="${await input.inputValue()}"`);

  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(300);
  const running = await page.evaluate(() => ({
    pane: !!document.querySelector('.u2-designer-pane'),
    running: document.querySelector('.u2-designer').classList.contains('u2-designer-running'),
    adorner: getComputedStyle(document.querySelector('.u2-designer-selected')).display,
  }));
  ok('view-and-panel/5b/run-mode-drops-glass-pane', !running.pane && running.running && running.adorner === 'none',
    JSON.stringify(running));

  await input.click();
  await input.press('End');
  await page.keyboard.type('ZZZ', {delay: 20});
  const typed = await input.inputValue();
  ok('view-and-panel/5c/run-mode-accepts-typing', typed !== before && /ZZZ/.test(typed), `"${before}" -> "${typed}"`);

  await clearBalloons(page);
  await page.locator('.u2-designer-surface [data-u2-name="saveButton"]').first().click();
  await page.waitForFunction(() => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => /cmd:save ran/.test(b.textContent)), null, {timeout: 5000}).catch(() => {});
  const fired = await balloons(page);
  await shot(page, 'view-and-panel-5-run-mode');
  ok('view-and-panel/5d/command-fires-in-run-mode', fired.some((b) => /cmd:save ran/.test(b)), JSON.stringify(fired));

  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
  ok('view-and-panel/5e/back-to-design', await page.evaluate(() => !!document.querySelector('.u2-designer-pane')));
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
  await shot(page, 'view-and-panel-6-reclick');
  ok('view-and-panel/6a/same-row-reclick-recovers', back === 'notesInput (u2-text-area)', `caption="${back}"`);
  ok('view-and-panel/6b/status-agrees', /notesInput/.test(await statusText(page)), await statusText(page));
}

async function checkOpen(page) {
  const prefilled = await openSpec(page, ROUND_TRIP);
  ok('view-and-panel/7a/open-dialog-prefilled',
    prefilled.includes('"$schema": "dg-ui/1"') && prefilled.includes('layout'),
    `${prefilled.length} chars`);
  // the selection is read before the tree is expanded: expanding a row selects it
  const status = await waitStatus(page, 'roundTrip');
  await expandTree(page);
  await shot(page, 'view-and-panel-7-round-trip');
  const after = await page.evaluate(() => ({
    named: [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name),
    value: document.querySelector('.u2-designer-surface [data-u2-name="rtName"] input')?.value,
    labels: [...document.querySelectorAll('.u2-designer-toolbox .u2-tree-label')].map((l) => l.textContent),
  }));
  ok('view-and-panel/7b/canvas-rerenders', after.value === 'Ethanol' && after.named.includes('rtFlag'),
    `value=${after.value} named=${JSON.stringify(after.named)}`);
  ok('view-and-panel/7c/tree-and-status-follow',
    after.labels.includes('rtForm') && /^roundTrip.* · 5 nodes/.test(status),
    `labels=${JSON.stringify(after.labels)} status="${status}"`);
  // the status bar is u2's authority on the selection (the platform panel renders one
  // current-object change behind and a lone programmatic set never renders — a known
  // platform defect; the click below is the user gesture that converges it, as elsewhere)
  const opened = await selectRow(page, 'roundTrip');
  ok('view-and-panel/7d/root-selected', /^roundTrip · /.test(status) && opened === 'roundTrip (u2-panel)',
    `status="${status}" caption-after-click="${opened}"`);

  await clearBalloons(page);
  await openSpec(page, '{"nope": 1}');
  await page.waitForTimeout(600);
  const refused = await balloons(page);
  const kept = await page.evaluate(() =>
    !!document.querySelector('.u2-designer-surface [data-u2-name="rtName"]'));
  ok('view-and-panel/7e/invalid-spec-refused', kept && refused.some((b) => /\$schema/.test(b)),
    `canvas kept=${kept} balloons=${JSON.stringify(refused)}`);
}

export const checks = [
  {id: 'view-and-panel/1 view renders', run: checkView},
  {id: 'view-and-panel/2 selection to context panel', run: checkPanel},
  {id: 'view-and-panel/3 canvas hit-test', run: checkHitTest},
  {id: 'view-and-panel/4 adorner geometry', run: checkAdorner},
  {id: 'view-and-panel/5 design inertness and run liveness', run: checkModes},
  {id: 'view-and-panel/6 same-row re-click', run: checkReclick},
  {id: 'view-and-panel/7 open round-trip', run: checkOpen},
];
