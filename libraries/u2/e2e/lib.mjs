/* What every designer check file shares: the readers (status bar, caption, panel, tree, canvas,
   chips), the drivers (select, drag, open, menus, pickers, the platform's own drag) and the two
   specs more than one file opens. Nodes are addressed by `data-u2-name` (WO-T1) wherever the spec
   names them; the platform chrome is addressed by its own classes. */
import {openApp} from './local.mjs';

/** What U2Demo's `u2DesignerApp` opens on (packages/U2Demo/src/designer.ts). */
export const APP = {package: 'U2demo', func: 'u2DesignerApp'};

const ROW = '.u2-designer-toolbox .u2-tree-row';
export const CAPTION = '.u2-designer-caption';
/** The status bar panel carries no class of its own — it is the leaf div reporting the node count. */
export const statusText = (page) => page.evaluate(() => [...document.querySelectorAll('div')]
  .find((x) => x.children.length === 0 && /\d+ node/.test(x.textContent))?.textContent ?? '');

export const caption = (page) => page.evaluate((sel) =>
  document.querySelector(sel)?.textContent ?? '', CAPTION);

export const balloons = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .map((b) => b.textContent.trim()).filter((t) => t));

export const clearBalloons = (page) => page.evaluate(() =>
  document.querySelectorAll('.d4-balloon-content, .d4-balloon').forEach((b) => b.remove()));

/** The property panel as `{section: {row: value}}` — what the context panel actually shows, read
 * off the read-only tables and the per-section form fields (stamped `data-u2-prop`) alike. */
export const panel = (page) => page.evaluate(() => {
  const props = document.querySelector('.u2-designer-properties');
  const sections = {};
  let current = null;
  const shown = (field) => {
    const input = field.querySelector('input, textarea, select');
    if (!input)
      return field.textContent.trim();
    return input.type === 'checkbox' ? String(input.checked) : input.value;
  };
  const scan = (el, into) => {
    for (const tr of el.querySelectorAll('tr'))
      into[tr.children[0]?.textContent.trim()] = tr.children[1]?.textContent.trim();
    for (const field of el.querySelectorAll('[data-u2-prop]'))
      into[field.dataset.u2Prop] = shown(field);
  };
  for (const el of props?.children ?? []) {
    if (el.tagName === 'H3')
      sections[current = el.textContent] = {};
    // the wiring sections are accordion panes — folded, their fields are still in the DOM
    else if (el.classList.contains('u2-accordion')) {
      current = null;
      for (const pane of el.querySelectorAll('.u2-accordion-pane'))
        scan(pane, sections[pane.querySelector('.u2-accordion-title')?.textContent ?? ''] = {});
    } else if (current)
      scan(el, sections[current]);
  }
  return sections;
});

/** Opens a folded panel section — Bindings starts closed on a node that binds nothing, so this is
 * the click a user makes before reaching for the `…` beside a path. */
export async function openSection(page, title) {
  const header = page.locator('.u2-designer-properties .u2-accordion-header')
    .filter({has: page.locator('.u2-accordion-title', {hasText: new RegExp(`^${title}$`)})}).first();
  if (await header.getAttribute('aria-expanded') === 'false')
    await header.click();
  await page.waitForTimeout(100);
}

/** The editor of one panel field, wherever it sits — a read-only field is disabled, so a name the
 * Properties and Bindings sections share resolves to the editable one. */
export const panelField = (page, name) => page.locator(`.u2-designer-properties [data-u2-prop="${name}"]`)
  .locator('input:not([disabled]), textarea:not([disabled]), select:not([disabled])').first();

/** The panel is asynchronous — the shell renders it after `grok.shell.o` is set. */
export async function waitCaption(page, name, timeout = 8000) {
  await page.waitForFunction(({sel, name}) =>
    (document.querySelector(sel)?.textContent ?? '').startsWith(`${name} (`),
  {sel: CAPTION, name}, {timeout}).catch(() => {});
  return caption(page);
}

/** The designer's own status bar: written synchronously by the selection, so it is the authority
 * on what the designer selected — the panel is the platform rendering that, one step behind. */
export async function waitStatus(page, name, timeout = 3000) {
  await page.waitForFunction((name) => [...document.querySelectorAll('div')]
    .some((x) => x.children.length === 0 && new RegExp(`${name} · \\d+ node`).test(x.textContent)),
  name, {timeout}).catch(() => {});
  return statusText(page);
}

/** One twistie at a time: expanding re-renders the virtual list, which invalidates every other
 * handle held across the click. */
export async function expandTree(page) {
  const collapsed = page.locator('.u2-designer-toolbox .u2-tree-twistie:not(.u2-tree-twistie-hidden)')
    .filter({hasText: '▸'});
  for (let i = 0; i < 30 && await collapsed.count() > 0; i++) {
    await collapsed.first().click({timeout: 2000}).catch(() => {});
    await page.waitForTimeout(80);
  }
}

export const row = (page, name) =>
  page.locator(ROW).filter({has: page.locator('.u2-tree-label', {hasText: name})}).first();

/** Selects a node from the structure tree. The platform's context panel runs one `grok.shell.o`
 * assignment behind (verified: reading `grok.shell.o` right after a single write still answers the
 * previous object), so a selection that has reached the status bar but not the panel is re-asserted
 * with the click a user would make — view-and-panel/6 covers that path on its own. */
export async function selectRow(page, name) {
  await row(page, name).click();
  await waitStatus(page, name);
  const shown = await waitCaption(page, name, 1500);
  if (shown.startsWith(`${name} (`))
    return shown;
  await row(page, name).click();
  return waitCaption(page, name);
}

export const ROUND_TRIP = JSON.stringify({
  $schema: 'dg-ui/1',
  root: {tag: 'u2-panel', name: 'roundTrip', children: [
    {tag: 'h2', props: {text: 'Round trip'}},
    {tag: 'u2-form', name: 'rtForm', children: [
      {tag: 'u2-text-input', name: 'rtName', props: {label: 'Reagent', value: 'Ethanol'}},
      {tag: 'u2-bool-input', name: 'rtFlag', props: {label: 'Checked', value: false}},
    ]},
  ]},
}, null, 2);

/** A property row, a per-child prop, a plain HTML node and a binding — `$.reagent` and `cmd:save`
 * resolve against the demo context the designer opened with. */
export const EDIT_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  root: {tag: 'u2-panel', name: 'editRoot', children: [
    {tag: 'h2', name: 'heading', props: {text: 'Editing'}},
    {tag: 'u2-tabs', name: 'tabs', children: [
      {tag: 'u2-panel', name: 'firstPane', props: {title: 'First'}, children: [
        {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}, bind: {value: '$.reagent'}},
      ]},
      {tag: 'u2-panel', name: 'secondPane', props: {title: 'Second'}},
    ]},
    {tag: 'u2-button', name: 'saveButton', props: {text: 'Save'}, on: {click: 'cmd:save'}},
  ]},
}, null, 2);

/** The Open button is one menu (Samples / Gallery / Open spec… / Save to gallery…) — this walks
 * it to the paste dialog. */
async function openSpecDialog(page) {
  await page.locator('.d4-ribbon').getByRole('button', {name: 'Open', exact: true}).first().click();
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu .u2-menu-label').filter({hasText: /^Open spec…$/}).first().click();
  await page.waitForSelector('.u2-designer-spec-editor textarea', {timeout: 10000});
}

/** Opens the dialog, answers what it was prefilled with — the designer's own dump — and commits
 * either the spec given or, with none, that same dump: the round trip. */
export async function openSpec(page, json) {
  await openSpecDialog(page);
  const area = page.locator('.u2-designer-spec-editor textarea');
  const prefilled = await area.inputValue();
  if (json !== undefined)
    await area.fill(json);
  await page.locator('.d4-dialog button, .u2-dialog button').filter({hasText: /^OK$/i}).first().click();
  return prefilled;
}

/** The Open dialog doubles as the dump reader: prefilled with the designer's own dump, cancelled
 * so the document — and its `modified` flag — stay untouched. */
export async function dumpViaDialog(page) {
  await openSpecDialog(page);
  const value = await page.locator('.u2-designer-spec-editor textarea').inputValue();
  await page.locator('.d4-dialog button, .u2-dialog button').filter({hasText: /^cancel$/i}).first().click();
  await page.waitForTimeout(200);
  return value;
}

export async function confirmDiscard(page) {
  await page.waitForSelector('.u2-dialog', {timeout: 5000});
  const title = await page.evaluate(() =>
    document.querySelector('.u2-dialog .u2-dialog-title-text')?.textContent ?? '');
  await page.locator('.u2-dialog button').filter({hasText: /^OK$/i}).first().click();
  await page.waitForTimeout(300);
  return title;
}

export const named = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name));

export const labels = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.u2-designer-toolbox .u2-tree-label')].map((l) => l.textContent));

export const nodeCount = async (page) => Number((await statusText(page)).match(/(\d+) nodes?/)?.[1] ?? -1);

export const ribbon = (page, name) => page.locator('.d4-ribbon').getByRole('button', {name, exact: true}).first();

export const paletteItem = (page, tag) => page.locator(`.u2-palette-item[data-u2-tag="${tag}"]`).first();

export const surfaceNode = (page, name) => page.locator(`.u2-designer-surface [data-u2-name="${name}"]`).first();

/** Many steps, not one jump: the drag layer starts on the 4px threshold and reads the drop target
 * off the mousemoves in between — a single move would land on neither. */
export async function drag(page, from, to, dy = 0) {
  const a = await from.boundingBox();
  const b = await to.boundingBox();
  await page.mouse.move(a.x + a.width / 2, a.y + a.height / 2);
  await page.mouse.down();
  await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2 + dy, {steps: 25});
  await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2 + dy + 1, {steps: 3});
  await page.mouse.up();
  await page.waitForTimeout(300);
}

/** The keyboard shortcuts are the view's own: they listen on the designer root and the toolbox, so
 * the focus has to be inside one of them — a click on a tree row is how a user puts it there. */
export async function focusTree(page, name = 'roundTrip') {
  await row(page, name).click();
  await page.waitForTimeout(100);
}

export async function toMode(page, mode) {
  await page.locator('.d4-ribbon').getByText(mode, {exact: true}).first().click();
  await page.waitForTimeout(500);
}

/** Starts over. Earlier checks leave the document modified, and `New form` confirms before it
 * discards — a clean document skips the dialog, so it is waited for, not required. */
export async function newForm(page) {
  await ribbon(page, 'New form').click();
  const dialog = page.locator('.u2-dialog').filter({hasText: 'Discard changes?'});
  await dialog.first().waitFor({state: 'visible', timeout: 2000}).catch(() => {});
  if (await dialog.count() > 0)
    await dialog.first().locator('button').filter({hasText: /^OK$/i}).first().click();
  await waitStatus(page, 'form1');
  await page.waitForTimeout(300);
}

/** Drops a palette item on a canvas node. The palette is filtered first: an unfiltered item sits
 * below the fold, and a drag that starts off-screen lands nowhere. */
export async function dropControl(page, tag, filter, onto) {
  await page.locator('.u2-palette input').first().fill(filter);
  await page.waitForTimeout(200);
  await drag(page, paletteItem(page, tag), surfaceNode(page, onto));
  await page.waitForTimeout(200);
}

/** Closes every view and applies the app function again — a second `u2DesignerApp()` over the very
 * object literal the first one was given. */
export async function reopenApp(page) {
  await page.evaluate(() => {
    for (const view of [...grok.shell.views])
      view.close();
  });
  await page.waitForTimeout(1000);
  await openApp(page, APP.package, APP.func);
  await page.waitForSelector('.u2-designer', {timeout: 60000});
  await page.waitForTimeout(1500);
}

const chip = (page, name) => page.locator(`.u2-designer-chip[data-u2-name="${name}"]`).first();

export const chipNames = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.u2-designer-chip')].map((c) => c.dataset.u2Name));

export const specErrors = (page) => page.evaluate(() =>
  document.querySelectorAll('.u2-designer-surface .u2-spec-error').length);

/** A chip is what a tray component is selected by — the same one-behind panel recovery a tree row
 * needs (selectRow), since the platform renders `grok.shell.o` a step late. */
export async function selectChip(page, name) {
  await chip(page, name).click();
  await waitStatus(page, name);
  const shown = await waitCaption(page, name, 1500);
  if (shown.startsWith(`${name} (`))
    return shown;
  await chip(page, name).click();
  return waitCaption(page, name);
}

/** The tray chip's own context menu: Refresh is the one design-time run a non-query ever gets. */
export async function chipMenu(page, name, label) {
  await chip(page, name).click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu .u2-menu-label').filter({hasText: new RegExp(`^${label}$`)}).first().click();
  await page.waitForTimeout(600);
}

/** Commits a guarded panel field (bindings/events commit on change, not per keystroke). */
export async function setPanelField(page, name, value) {
  const field = panelField(page, name);
  await field.fill(value);
  await field.press('Enter');
  await page.waitForTimeout(400);
}

export const pickerRow = (page, label) => page.locator('.u2-bind-picker .u2-list-row')
  .filter({has: page.locator('.u2-tree-label', {hasText: label})}).first();

export const pickerLabels = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.u2-bind-picker .u2-tree-label')].map((l) => l.textContent));

/** The path field of one binding row. A prop name lives in the Properties section too, and an
 * unbound prop's editor there is enabled — the picker button beside it is what tells them apart. */
export const bindField = (page, name) => page.locator(`.u2-designer-properties [data-u2-prop="${name}"]`)
  .filter({has: page.locator(`[data-u2-bind-pick="${name}"]`)}).locator('input').first();

const fieldOf = (page, name) =>
  page.locator(`.u2-designer-surface [data-u2-name="${name}"] input`).first();

export const fieldValue = (page, name) => fieldOf(page, name).inputValue();

export async function setField(page, name, value) {
  await fieldOf(page, name).fill(value);
  await page.waitForTimeout(400);
}

export const trayMenu = async (page, ...labels) => {
  await page.locator('.u2-designer-tray-add').first().click();
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  for (const label of labels) {
    await page.locator('.u2-menu .u2-menu-label').filter({hasText: new RegExp(`^${label}$`)}).first().click();
    await page.waitForTimeout(200);
  }
  await page.waitForTimeout(400);
};

const funcRow = (page, name) => page.locator(`.u2-func-picker .u2-list-row[data-u2-func="${name}"]`).first();

export const dialogOK = (page, dialog) =>
  page.locator(`${dialog} button`).filter({hasText: /^OK$/}).first().click();

export const dialogCancel = (page, dialog) =>
  page.locator(`${dialog} button`).filter({hasText: /^CANCEL$/}).first().click();

export async function pickFunc(page, name) {
  await page.waitForSelector('.u2-func-picker', {timeout: 5000});
  await page.locator('.u2-func-picker input').first().fill(name);
  await page.waitForTimeout(300);
  await funcRow(page, name).click();
  await page.waitForTimeout(300);
}

/** Expands the picker down a chain of labels; answers what the last expansion showed. */
export async function pickerWalk(page, ...labels) {
  await page.waitForSelector('.u2-bind-picker', {timeout: 5000});
  for (const label of labels) {
    await pickerRow(page, label).locator('.u2-tree-twistie').click();
    await page.waitForTimeout(300);
  }
  return pickerLabels(page);
}

/** Binds through the picker: walks the chain, then picks the row whose label starts with `leaf`.
 * Answers what the tree showed, so a failure names it. */
export async function pickLeaf(page, leaf, ...walk) {
  const labels = await pickerWalk(page, ...walk);
  const found = labels.find((label) => label.startsWith(leaf));
  if (found !== undefined) {
    await pickerRow(page, found).click();
    await dialogOK(page, '.u2-bind-picker');
    await page.waitForTimeout(400);
  }
  return {found, labels};
}

/** Selects a control and binds one of its props the way a user does: the `…` beside the Bindings
 * row, then the picker walked down to the leaf. */
export async function bindThroughPicker(page, control, leaf, ...walk) {
  await selectRow(page, control);
  await openSection(page, 'Bindings');
  await page.locator('.u2-designer-properties [data-u2-bind-pick="value"]').first().click();
  return pickLeaf(page, leaf, ...walk);
}

/** A platform function whose firing is observable from the page. U2Demo publishes `u2Record` — it
 * takes an argument and leaves a durable trace (the designer's Run log; the balloon is gone in five
 * seconds) — so the call itself is what a check reads. Where the local fixture has not caught up
 * with it, the session registers a stub of the same signature instead. */
export async function ensureRecorder(page) {
  await page.evaluate(() => {
    window.__u2Fired = '';
    if (DG.Func.find({name: 'u2Record'}).length === 0) {
      grok.functions.register({signature: 'String u2Record(string text)',
        run: (text) => window.__u2Fired = text, isAsync: false});
      return;
    }
    window.__u2RecordSpy?.unsubscribe?.();
    window.__u2RecordSpy = grok.functions.onBeforeRunAction.subscribe((call) => {
      if ((call?.func?.name ?? '') === 'u2Record')
        window.__u2Fired = call?.inputs?.text ?? '';
    });
  });
}

/** A table the user already has open — what the workspace picker and a `u2-table-source` read. */
export const ensureDemog = (page) => page.evaluate(() => {
  if (!grok.shell.tableNames.includes('u2demog')) {
    const table = grok.data.demo.demog(50);
    table.name = 'u2demog';
    // a table nothing has a view on has no current row, and a bound field would show nothing
    table.currentRowIdx = 0;
    grok.shell.addTable(table);
  }
});

/** Every function the platform runs from here on — the design-time silence proof is a run spy,
 * not a network spy: in local mode a client-side function makes no request either way. */
export const spyRuns = (page) => page.evaluate(() => {
  window.__u2Sub?.unsubscribe?.();
  window.__u2Runs = [];
  window.__u2Sub = grok.functions.onBeforeRunAction
    .subscribe((call) => window.__u2Runs.push(call?.func?.name ?? '?'));
});

/** Arms a drag source carrying `spec` and drives the gesture onto `target`. `mid` runs with the
 * button still down, past the 5px threshold. What the source carries:
 * `{path}` a FileInfo built from the string, `{dapiPath}` the server's own FileInfo for that path,
 * `{paths}` an array of them (what a card-view multi-selection drags), `{func}` a registered func. */
export async function platformDrag(page, spec, target, mid) {
  await page.evaluate(async (spec) => {
    document.querySelector('#u2-drag-source')?.remove();
    const el = document.createElement('div');
    el.id = 'u2-drag-source';
    el.textContent = 'drag';
    el.style.cssText = 'position:fixed;left:8px;top:430px;width:90px;height:26px;' +
      'z-index:100000;background:#ddd';
    document.body.append(el);
    const real = async (full) => (await grok.dapi.files.list(full.replace(/\/[^/]*$/, '/')))
      .find((f) => f.fullPath === full);
    let object = DG.Func.find({name: spec.func})[0];
    if (spec.dapiPath)
      object = await real(spec.dapiPath);
    else if (spec.paths)
      object = await Promise.all(spec.paths.map(real));
    else if (spec.path)
      object = DG.FileInfo.fromString(spec.path, '');
    if (object === undefined || (Array.isArray(object) && object.some((x) => x === undefined)))
      throw new Error(`drag source resolved to nothing: ${JSON.stringify(spec)}`);
    ui.makeDraggable(el, {getDragObject: () => object});
  }, spec);
  const from = await page.locator('#u2-drag-source').boundingBox();
  const to = await target.boundingBox();
  await page.mouse.move(from.x + from.width / 2, from.y + from.height / 2);
  await page.mouse.down();
  await page.mouse.move(to.x + to.width / 2, to.y + to.height / 2, {steps: 25});
  await page.mouse.move(to.x + to.width / 2, to.y + to.height / 2 + 1, {steps: 3});
  if (mid)
    await mid();
  await page.mouse.up();
  await page.waitForTimeout(800);
  await page.evaluate(() => document.querySelector('#u2-drag-source')?.remove());
}

/** What the platform leaves on screen for the duration of a drag, and must take away after it. */
export const dragChrome = (page) => page.evaluate(() => ({
  body: document.body.classList.contains('d4-drag'),
  ghosts: document.querySelectorAll('.d4-column-drag-host').length,
  tray: document.querySelector('.u2-designer-tray')?.classList.contains('u2-designer-tray-drop'),
  hover: getComputedStyle(document.querySelector('.u2-designer-hover')).display,
}));

export const components = async (page) => JSON.parse(await dumpViaDialog(page)).components ?? [];
