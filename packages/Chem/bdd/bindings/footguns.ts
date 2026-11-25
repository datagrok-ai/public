// Foot-gun flags (features/footguns: the scenarios of what each check expects): the middle slot of an editor
// dialog's footer (src/utils/footguns/footer-slot.ts) and the context panel's Checks section
// (src/widgets/checks-panel.ts), with what the user does to raise them: a value put into the Crux sketcher of the cell
// editor, a reaction typed into New Reaction, a reaction picked in a Transform dialog, a table of molecules or
// templates. The slot is `.chem-footgun-slot` (its level in `data-level`, its whole text in `aria-label`); the Checks
// section's lines are `.chem-footgun-checks > [data-level]`. Each step fails on a build without them: a step that
// finds no slot fails, and "shows no flag" asks for a slot that is there and empty in Crux mode.
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, expect, pickTopMenu, pollMs} from '@datagrok-libraries/bdd/runtime';
import {clickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';

declare const grok: any;
declare const DG: any;

type Rect = {left: number, top: number, width: number, height: number};

/** The reaction cell editor's dialog (src/utils/reactions/reaction-cell-editor.ts): the one holding its box. */
const REACTION_EDITOR = '.d4-dialog:has(.chem-reaction-cell-editor)';

/** The dialog a step names, as a CSS selector: the cell editor's (the one with a molecule sketcher in it), one by its
 * title, or the topmost. */
function dialogSelector(which: string | null): string {
  if (which === 'sketcher')
    return '.d4-dialog:has(.chem-sketcher-host, input[placeholder^="SMILES"])';
  if (which === 'reaction')
    return REACTION_EDITOR;
  return '.d4-dialog';
}

/** In the page: the dialog a step names (the last match: the topmost), or null. */
async function hasDialog(page: Page, which: string | null, title?: string): Promise<boolean> {
  return page.evaluate(([sel, t]) => {
    const all = [...document.querySelectorAll(sel)].filter((d) => !t ||
      (d.querySelector('.d4-dialog-title, .u2-dialog-title-text')?.textContent ?? '').trim() === t);
    return all.length > 0;
  }, [dialogSelector(which), title ?? ''] as [string, string]);
}

/** In the page: `fn` run on the dialog a step names, once it is open (the editors open asynchronously: the cell editor
 * after a double-click, New Reaction after its button); it throws when none opens within 30 s. */
async function inDialog<T>(page: Page, which: string | null, title: string | undefined, fn: string): Promise<T> {
  await expect.poll(() => hasDialog(page, which, title), {message: `the ${title ? `"${title}"` : which ?? ''} dialog open`,
    timeout: pollMs(30000)}).toBe(true);
  return page.evaluate(([sel, t, body]) => {
    const all = [...document.querySelectorAll(sel)].filter((d) => !t ||
      (d.querySelector('.d4-dialog-title, .u2-dialog-title-text')?.textContent ?? '').trim() === t);
    const dialog = all[all.length - 1] as HTMLElement | undefined;
    if (!dialog)
      throw new Error(`no ${t ? `"${t}"` : sel} dialog is open`);
    return new Function('dialog', body)(dialog);
  }, [dialogSelector(which), title ?? '', fn] as [string, string, string]) as Promise<T>;
}

/** The dialog's footer buttons by their text (OK, SAVE, CANCEL), on screen. */
async function buttonRects(page: Page, which: string | null, title: string | undefined, names: string[]):
  Promise<{[name: string]: Rect}> {
  return inDialog(page, which, title, `
    const out = {};
    for (const name of ${JSON.stringify(names)}) {
      const b = [...dialog.querySelectorAll('.d4-dialog-footer button')]
        .find((x) => x.textContent.trim().toUpperCase() === name.toUpperCase());
      if (!b)
        throw new Error('no ' + name + ' button in the dialog footer');
      const r = b.getBoundingClientRect();
      out[name] = {left: r.left, top: r.top, width: r.width, height: r.height};
    }
    return out;`);
}

/** What the footer's middle slot shows: whether it is there, its level, the text it shows and the whole text, its box,
 * and whether its text is cut. */
interface SlotState {present: boolean, level: string, text: string, whole: string, rect: Rect | null, cut: boolean}

async function slotState(page: Page, which: string | null, title?: string): Promise<SlotState> {
  return inDialog(page, which, title, `
    const slot = dialog.querySelector('.d4-dialog-footer .chem-footgun-slot');
    if (!slot)
      return {present: false, level: '', text: '', whole: '', rect: null, cut: false};
    const label = slot.lastElementChild;
    const r = slot.getBoundingClientRect();
    return {present: true, level: slot.dataset.level ?? '', text: slot.textContent.trim(),
      whole: slot.getAttribute('aria-label') ?? '', rect: {left: r.left, top: r.top, width: r.width, height: r.height},
      cut: !!label && label.scrollWidth > label.clientWidth + 1};`);
}

/** The page's engine (its Chem setting as Chem reads it in this page). */
async function engineOf(page: Page): Promise<string> {
  return page.evaluate(() => grok.functions.call('Chem:chemicalEngine'));
}

/** Waits for the packages' autostarts (Chem's menus and panels are there only after them: before, a step could find
 * "no Chem in the top menu"). */
async function autostarted(page: Page): Promise<void> {
  await page.evaluate(async () => {
    await grok.shell.autostartsCompleted;
  });
}

/** The noted button positions, per page. */
const noted = new WeakMap<Page, {[name: string]: Rect}>();

/** The values each column of the last table a step opened held as it opened, per page (for "exactly what it was"). */
const opened = new WeakMap<Page, {[column: string]: string[]}>();

// ---------------------------------------------------------------------------------------------------- tables

export const moleculeTable = Given('user opens a table {string} with a {word} column {string} holding {string}',
  async (page: Page, name: string, kind: string, column: string, value: string) => {
    await autostarted(page);
    await page.evaluate(async ([n, k, c, v]) => {
      const molblock = k.toLowerCase() === 'molblock';
      const cell = molblock ? grok.chem.convert(v, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock) : v;
      const col = DG.Column.fromStrings(c, [cell]);
      col.semType = DG.SEMTYPE.MOLECULE;
      col.meta.units = molblock ? DG.chem.Notation.MolBlock : DG.chem.Notation.Smiles;
      const df = DG.DataFrame.fromColumns([col]);
      df.name = n;
      grok.shell.addTableView(df);
    }, [name, kind, column, value] as [string, string, string, string]);
    await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid &&
      !grok.shell.tv.grid.isRenderPending, name), {message: `a table view of "${name}"`}).toBe(true);
    // the value as the table opened, for "exactly what it was" (the engine-poison scenarios)
    opened.set(page, {[column]: [await page.evaluate((c) => String(grok.shell.tv.dataFrame.col(c).get(0)), column)]});
  }, {tier: 'api', description: 'a one-row table of a molecule column of the units named (SMILES, or molblock: the SMILES converted to a V2000 molblock), in a table view of its own'});

/** Opens a table of one reaction column holding `values`, in a table view of its own, and notes the values (for
 * "exactly what it was"). */
async function openReactionTable(page: Page, name: string, column: string, units: string, values: string[]):
  Promise<void> {
  await autostarted(page);
  await page.evaluate(([n, c, u, v]) => {
    const col = DG.Column.fromStrings(c, v);
    col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
    col.meta.units = u;
    col.meta.cellRenderer = 'ChemicalReaction';
    const df = DG.DataFrame.fromColumns([col]);
    df.name = n;
    grok.shell.addTableView(df);
  }, [name, column, units, values] as [string, string, string, string[]]);
  opened.set(page, {[column]: values});
  await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid &&
    !grok.shell.tv.grid.isRenderPending, name), {message: `a table view of "${name}"`}).toBe(true);
}

export const reactionTable = Given('user opens a table {string} with a reaction column {string} of units {string}:',
  async (page: Page, name: string, column: string, units: string, rows: string[][]) => {
    const [, ...body] = rows;
    await openReactionTable(page, name, column, units, body.map((r) => r[0]));
  }, {tier: 'api', description: 'a table of one reaction column (the header row names it; its values as written), its semantic type ChemicalReaction and its units as named ("smarts": a template column), in a table view of its own'});

export const detectedTable = Given('user opens a CSV {string} with a column {string} whose semantic types are detected:',
  async (page: Page, name: string, column: string, rows: string[][]) => {
    const [, ...body] = rows;
    await autostarted(page);
    await page.evaluate(async ([n, c, values]) => {
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(c, values)]);
      df.name = n;
      await grok.data.detectSemanticTypes(df);
      grok.shell.addTableView(df);
    }, [name, column, body.map((r) => r[0])] as [string, string, string[]]);
    opened.set(page, {[column]: body.map((r) => r[0])});
    await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid &&
      !grok.shell.tv.grid.isRenderPending, name), {message: `a table view of "${name}"`}).toBe(true);
  }, {tier: 'api', description: 'a table of one text column (the header row names it; its values as written) whose semantic type and units the platform\'s detectors give it (Chem\'s reaction detector among them), as for an opened CSV, in a table view of its own'});

export const pageEngine = Given('the Chemical Engine of this page is {string}', async (page: Page, engine: string) => {
  // Chem reads its settings from its package (re-read at most every second): this page's reads give the engine named
  // until the shell reset that ends the scenario (the page runs the worker's next scenarios too: left in place, the
  // override made every later scenario's "the Chemical Engine is Crux" gate skip); the server's stored settings are
  // not touched
  await page.evaluate((e) => {
    const w = window as any;
    const proto = DG.Package.prototype;
    const original = w.__guardSettings ?? Object.getOwnPropertyDescriptor(proto, 'settings')!;
    w.__guardSettings = original;
    Object.defineProperty(proto, 'settings', {configurable: true, get(this: any) {
      const s = original.get!.call(this);
      return this.name === 'Chem' ? {...(s ?? {}), ChemicalEngine: e} : s;
    }});
    const shell = grok.shell;
    if (!w.__guardCloseAll) {
      w.__guardCloseAll = true;
      // the shell reset's closeAll puts the server's settings back first
      shell.closeAll = function(this: any, ...args: any[]) {
        Object.defineProperty(proto, 'settings', w.__guardSettings);
        delete w.__guardSettings;
        delete shell.closeAll;
        delete w.__guardCloseAll;
        // Chem reads its settings at most once a second: the next scenario's first read must be the server's
        const until = performance.now() + 1100;
        while (performance.now() < until) { /* wait out Chem's settings cache */ }
        return shell.closeAll(...args);
      };
    }
  }, engine);
  await expect.poll(() => engineOf(page), {message: `Chem in this page runs ${engine}`, timeout: pollMs(15000)})
    .toBe(engine);
}, {tier: 'api', description: 'Chem in this page reads its Chemical Engine setting as the engine named (and redraws as at a switch); the stand\'s stored settings are untouched'});

// ---------------------------------------------------------------------------------------------------- the sketcher

/** Puts `value` into the Crux sketcher of the topmost sketcher dialog through its host (the platform's Sketcher), as a
 * SMILES or as a query (SMARTS), once Crux is ready; Crux then holds it as if drawn, and the host as the value OK
 * writes. */
async function putIntoCrux(page: Page, value: string, asQuery: boolean): Promise<void> {
  await expect.poll(() => page.evaluate(() => {
    const dialogs = [...document.querySelectorAll('.d4-dialog')].filter((d) => d.querySelector('.crux-sketcher'));
    const root = dialogs[dialogs.length - 1]?.querySelector('.crux-sketcher');
    const impl = root ? DG.Widget.find(root) : null;
    return !!impl?.isInitialized && !impl.isRenderPending;
  }), {message: 'Crux ready in the sketcher dialog', timeout: pollMs(60000)}).toBe(true);
  await page.evaluate(([v, q]) => {
    const dialogs = [...document.querySelectorAll('.d4-dialog')].filter((d) => d.querySelector('.crux-sketcher'));
    const root = dialogs[dialogs.length - 1]!.querySelector('.crux-sketcher')!;
    let host: any = null;
    for (let e = root.parentElement; e !== null && host === null; e = e.parentElement) {
      const w = DG.Widget.find(e);
      if (w && w.onSketcherReady)
        host = w;
    }
    if (!host)
      throw new Error('no sketcher host round Crux');
    q ? host.setSmarts(v) : host.setSmiles(v);
  }, [value, asQuery] as [string, boolean]);
  // Crux shows it (a value Crux can't read changes nothing: the step says so)
  await expect.poll(() => page.evaluate(() => {
    const dialogs = [...document.querySelectorAll('.d4-dialog')].filter((d) => d.querySelector('.crux-sketcher'));
    const crux: any = dialogs[dialogs.length - 1]?.querySelector('crux-sketch');
    return crux && !crux.isEmpty;
  }), {message: `Crux shows ${value}`, timeout: pollMs(15000)}).toBe(true);
}

export const putMolecule = When('user puts {string} into the Crux sketcher of the dialog', (page: Page, value: string) =>
  putIntoCrux(page, value, false),
{tier: 'api', description: 'the molecule set on the sketcher dialog\'s host as SMILES once Crux is ready: Crux shows it, and the host holds it as the value OK writes, as after drawing it'});

export const putQuery = When('user puts the query {string} into the Crux sketcher of the dialog',
  (page: Page, value: string) => putIntoCrux(page, value, true),
{tier: 'api', description: 'the query set on the sketcher dialog\'s host as SMARTS once Crux is ready: Crux shows it as a query, as after drawing an atom list or a query bond'});

export const resizeSketcherDialog = When('user resizes the sketcher dialog to {int} px wide',
  async (page: Page, width: number) => {
    const box = await inDialog<Rect>(page, 'sketcher', undefined, `
      const r = dialog.getBoundingClientRect();
      return {left: r.left, top: r.top, width: r.width, height: r.height};`);
    // the dialog's right edge, dragged
    const x = box.left + box.width - 2;
    const y = box.top + box.height / 2;
    await page.mouse.move(x, y);
    await page.mouse.down();
    await page.mouse.move(x - (box.width - width), y, {steps: 10});
    await page.mouse.up();
    await expect.poll(() => inDialog<number>(page, 'sketcher', undefined,
      'return Math.round(dialog.getBoundingClientRect().width);'), {message: `the sketcher dialog ${width} px wide`})
      .toBeLessThanOrEqual(width + 2);
  }, {tier: 'ui', description: 'the sketcher dialog\'s right edge dragged until the dialog is that wide'});

// ---------------------------------------------------------------------------------------------------- dialogs

export const clickDialogButton = When('user clicks the {string} button of the dialog', async (page: Page, text: string) => {
  const sel = '.d4-dialog';
  const button = page.locator(sel).last().locator('button', {hasText: text}).first();
  await expect(button, `the "${text}" button of the dialog`).toBeVisible({timeout: pollMs(30000)});
  await button.click();
}, {tier: 'ui', description: 'a click on the topmost dialog\'s button whose text holds the words given ("New Reaction" finds "+ New Reaction")'});

/** The New Reaction editor's input by its caption. */
async function fillEditorInput(page: Page, caption: string, value: string): Promise<void> {
  await expect.poll(() => hasDialog(page, null, 'New Reaction'), {message: 'the New Reaction dialog',
    timeout: pollMs(30000)}).toBe(true);
  const dialog = page.locator('.d4-dialog', {has: page.locator('.d4-dialog-title, .u2-dialog-title-text',
    {hasText: 'New Reaction'})}).last();
  const input = caption === 'Reaction SMARTS' ? dialog.locator('textarea').first() :
    dialog.locator('.ui-input-root', {has: page.locator('label', {hasText: caption})}).locator('input').first();
  await input.fill(value);
  await input.dispatchEvent('input');
}

export const enterSmarts = When('user enters {string} as the reaction SMARTS', (page: Page, smarts: string) =>
  fillEditorInput(page, 'Reaction SMARTS', smarts),
{tier: 'ui', description: 'the text typed into the New Reaction dialog\'s Reaction SMARTS box (it validates and draws as the user types)'});

export const enterName = When('user enters {string} as the reaction name', (page: Page, name: string) =>
  fillEditorInput(page, 'Reaction Name', name),
{tier: 'ui', description: 'the text typed into the New Reaction dialog\'s Reaction Name field'});

const noteButtons = async (page: Page, a: string, b: string, which: string | null, title?: string) => {
  // settled: the same boxes twice, 300 ms apart (Crux sizes its dialog as it loads)
  let last = '';
  await expect.poll(async () => {
    const now = JSON.stringify(await buttonRects(page, which, title, [a, b]));
    const same = now === last;
    last = now;
    return same;
  }, {message: 'the dialog\'s buttons settled', intervals: [300], timeout: pollMs(30000)}).toBe(true);
  noted.set(page, JSON.parse(last));
};

export const noteSketcherButtons = When('user notes where the {word} and {word} buttons of the sketcher dialog are',
  (page: Page, a: string, b: string) => noteButtons(page, a, b, 'sketcher'),
{tier: 'api', description: 'the boxes of the two footer buttons of the sketcher dialog, once they hold still, kept for a later step'});

export const noteNamedButtons = When('user notes where the {word} and {word} buttons of the {string} dialog are',
  (page: Page, a: string, b: string, title: string) => noteButtons(page, a, b, null, title),
{tier: 'api', description: 'the boxes of the two footer buttons of the dialog with that title, once they hold still, kept for a later step'});

const buttonsStill = async (page: Page, a: string, b: string, which: string | null, title?: string) => {
  const before = noted.get(page);
  if (!before)
    throw new Error('no button positions were noted earlier in the scenario');
  const now = await buttonRects(page, which, title, [a, b]);
  for (const name of [a, b]) {
    for (const k of ['left', 'top', 'width', 'height'] as const) {
      expect(Math.abs(now[name][k] - before[name][k]), `${name}'s ${k}: ${before[name][k]} px before, ${now[name][k]} now`)
        .toBeLessThanOrEqual(0.5);
    }
  }
};

export const sketcherButtonsStill = Then('the {word} and {word} buttons of the sketcher dialog are where they were',
  (page: Page, a: string, b: string) => buttonsStill(page, a, b, 'sketcher'),
{description: 'each button\'s box (left, top, width, height) within 0.5 px of the one noted'});

export const namedButtonsStill = Then('the {word} and {word} buttons of the {string} dialog are where they were',
  (page: Page, a: string, b: string, title: string) => buttonsStill(page, a, b, null, title),
{description: 'each button\'s box (left, top, width, height) within 0.5 px of the one noted'});

export const buttonsClearOfSlot = Then('the OK and CANCEL buttons do not overlap the middle slot', async (page: Page) => {
  const slot = await slotState(page, 'sketcher');
  expect(slot.present, 'a middle slot in the footer').toBe(true);
  const buttons = await buttonRects(page, 'sketcher', undefined, ['OK', 'CANCEL']);
  for (const [name, r] of Object.entries(buttons)) {
    const s = slot.rect!;
    const overlap = Math.min(r.left + r.width, s.left + s.width) - Math.max(r.left, s.left);
    expect(overlap, `${name} and the slot overlap by ${overlap} px`).toBeLessThanOrEqual(0);
  }
}, {description: 'the slot\'s box ends left of both buttons\' boxes'});

// ---------------------------------------------------------------------------------------------------- the slot

/** Polls the slot until `ok(state)`. */
async function expectSlot(page: Page, which: string | null, title: string | undefined, message: string,
  ok: (s: SlotState) => boolean): Promise<SlotState> {
  let state: SlotState | null = null;
  await withState(expect.poll(async () => {
    state = await slotState(page, which, title);
    return ok(state);
  }, {message, timeout: pollMs(15000)}).toBe(true), () => `${message}; the slot: ${JSON.stringify(state)}`);
  return state!;
}

/** Runs a poll; when it fails, its error starts with `said()` (what the poll last saw). */
async function withState(poll: Promise<unknown>, said: () => string): Promise<void> {
  try {
    await poll;
  } catch (e) {
    if (e instanceof Error)
      e.message = `${said()}\n${e.message}`;
    throw e;
  }
}

const hasWords = (text: string, words: string[]) => words.every((w) => text.toLowerCase().includes(w.toLowerCase()));

/** The slot's flags of `level` (its whole text, one line per flag). */
const linesOf = (s: SlotState) => s.whole.split('\n').filter((l) => l);

export const slotEmptySketcher = Then('the footer of the sketcher dialog has an empty middle slot', (page: Page) =>
  expectSlot(page, 'sketcher', undefined, 'an empty middle slot in the footer', (s) => s.present && s.text === ''),
{description: 'the footer holds its middle slot, and it says nothing'});

/** "shows no flag": in Crux mode the slot is there and empty; in RDKit mode there is none, or it is empty. */
const noFlag = async (page: Page, which: string | null, title?: string) => {
  const crux = await engineOf(page) === 'Crux';
  await expectSlot(page, which, title, crux ? 'the middle slot there and empty (Crux mode)' : 'no flag (RDKit mode)',
    (s) => (crux ? s.present : true) && s.text === '');
};

export const noFlagSketcher = Then('the footer of the sketcher dialog shows no flag', (page: Page) =>
  noFlag(page, 'sketcher'),
{description: 'Crux mode: the middle slot is there and empty; RDKit mode: no slot, or an empty one'});

export const noFlagNamed = Then('the footer of the {string} dialog shows no flag', (page: Page, title: string) =>
  noFlag(page, null, title),
{description: 'Crux mode: the middle slot is there and empty; RDKit mode: no slot, or an empty one'});

export const slotShows = Then('the middle slot of the footer shows a {word} flag', (page: Page, level: string) =>
  expectSlot(page, null, undefined, `a ${level} flag in the middle slot`,
    (s) => s.present && s.level === level && s.text !== ''),
{description: 'the topmost dialog\'s middle slot shows a flag, the worst of its flags of that level'});

export const slotShowsWord = Then('the middle slot of the footer shows a {word} flag mentioning {string}',
  (page: Page, level: string, word: string) => expectSlot(page, null, undefined,
    `a ${level} flag mentioning "${word}"`, (s) => s.present && s.level === level && hasWords(s.whole, [word])),
{description: 'the middle slot\'s worst flag is of that level, and its words hold the word (any case)'});

export const slotShowsWords = Then('the middle slot of the footer shows a {word} flag mentioning {string} and {string}',
  (page: Page, level: string, a: string, b: string) => expectSlot(page, null, undefined,
    `a ${level} flag mentioning "${a}" and "${b}"`, (s) => s.present && s.level === level &&
      linesOf(s).some((l) => hasWords(l, [a, b]))),
{description: 'the middle slot\'s worst flag is of that level, and one flag\'s words hold both words (any case)'});

export const slotNoBlock = Then('the middle slot of the footer shows no block flag', (page: Page) =>
  expectSlot(page, null, undefined, 'a middle slot without a block', (s) => s.present && s.level !== 'block'),
{description: 'the topmost dialog\'s middle slot is there, and shows no block'});

export const slotNoWord = Then('the middle slot of the footer shows no flag mentioning {string}',
  (page: Page, word: string) => expectSlot(page, null, undefined, `a middle slot with no flag mentioning "${word}"`,
    (s) => s.present && !linesOf(s).some((l) => hasWords(l, [word]))),
{description: 'the topmost dialog\'s middle slot is there (its checks ran as it opened); no flag of it holds the word'});

export const namedNoBlock = Then('the footer of the {string} dialog shows no block flag', (page: Page, title: string) =>
  expectSlot(page, null, title, 'a middle slot without a block', (s) => s.present && s.level !== 'block'),
{description: 'the dialog\'s middle slot is there, and shows no block'});

export const slotCut = Then('the middle slot of the footer shows a flag whose text is cut to the slot\'s width',
  async (page: Page) => {
    const s = await expectSlot(page, null, undefined, 'a flag in the middle slot', (x) => x.present && x.text !== '');
    expect(s.cut, 'the slot\'s text runs past the slot and is cut').toBe(true);
  }, {description: 'the slot shows a flag whose text is wider than the slot (cut with an ellipsis)'});

export const hoverSlot = When('user hovers over the middle slot of the footer', async (page: Page) => {
  const slot = page.locator('.d4-dialog').last().locator('.chem-footgun-slot');
  await expect(slot, 'the footer\'s middle slot').toBeVisible();
  await slot.hover();
}, {tier: 'ui', description: 'the pointer over the topmost dialog\'s middle slot'});

export const slotTooltip = Then('a tooltip shows the flag\'s whole text, naming every problem found', async (page: Page) => {
  const s = await slotState(page, null);
  const lines = linesOf(s);
  expect(lines.length, 'the slot holds more than one problem').toBeGreaterThan(1);
  const tooltip = page.locator('.d4-tooltip').last();
  await expect(tooltip, 'a tooltip').toBeVisible({timeout: pollMs(5000)});
  const text = (await tooltip.innerText()).replace(/\s+/g, ' ');
  for (const line of lines)
    expect(text, `the tooltip names "${line}"`).toContain(line.replace(/\s+/g, ' '));
}, {description: 'the tooltip holds every flag of the slot whole (more than the one it shows), one after another'});

// ---------------------------------------------------------------------------------------------------- OK / SAVE

const buttonDisabled = async (page: Page, name: string, which: string | null, title?: string) =>
  inDialog<boolean>(page, which, title, `
    const b = [...dialog.querySelectorAll('.d4-dialog-footer button')]
      .find((x) => x.textContent.trim().toUpperCase() === ${JSON.stringify(name.toUpperCase())});
    if (!b)
      throw new Error('no ${name} button');
    return b.disabled || b.classList.contains('d4-disabled') || b.classList.contains('disabled');`);

const sketcherButtonIs = (want: boolean) => async (page: Page, name: string) =>
  expect.poll(() => buttonDisabled(page, name, 'sketcher'), {message: `${name} ${want ? 'disabled' : 'enabled'}`})
    .toBe(want);
const namedButtonIs = (want: boolean) => async (page: Page, name: string, title: string) =>
  expect.poll(() => buttonDisabled(page, name, null, title), {message: `${name} ${want ? 'disabled' : 'enabled'}`})
    .toBe(want);
const topButtonIs = (want: boolean) => async (page: Page, name: string) =>
  expect.poll(() => buttonDisabled(page, name, null), {message: `${name} ${want ? 'disabled' : 'enabled'}`})
    .toBe(want);

export const sketcherButtonDisabled = Then('the {word} button of the sketcher dialog is disabled',
  sketcherButtonIs(true), {description: 'the sketcher dialog\'s footer button is disabled (a disabled OK neither clicks nor answers Enter)'});
export const sketcherButtonEnabled = Then('the {word} button of the sketcher dialog is enabled',
  sketcherButtonIs(false), {description: 'the sketcher dialog\'s footer button is enabled'});
export const namedButtonDisabled = Then('the {word} button of the {string} dialog is disabled',
  namedButtonIs(true), {description: 'the footer button of the dialog with that title is disabled'});
export const namedButtonEnabled = Then('the {word} button of the {string} dialog is enabled',
  namedButtonIs(false), {description: 'the footer button of the dialog with that title is enabled'});
export const topButtonDisabled = Then('the {word} button of the dialog is disabled',
  topButtonIs(true), {description: 'the topmost dialog\'s footer button is disabled'});
export const topButtonEnabled = Then('the {word} button of the dialog is enabled',
  topButtonIs(false), {description: 'the topmost dialog\'s footer button is enabled'});

export const clickNamedButton = When('user clicks on {word} button of the {string} dialog',
  async (page: Page, name: string, title: string) => {
    const dialog = page.locator('.d4-dialog', {has: page.locator('.d4-dialog-title, .u2-dialog-title-text',
      {hasText: title})}).last();
    await dialog.locator('.d4-dialog-footer button', {hasText: new RegExp(`^\\s*${name}\\s*$`, 'i')}).click();
  }, {tier: 'ui', description: 'a click on the footer button of the dialog with that title'});

// ---------------------------------------------------------------------------------------------------- cells

export const valueChanged = Then('the value in row {int} of {string} column is no longer {string}',
  async (page: Page, row: number, column: string, before: string) => {
    await expect.poll(() => page.evaluate(([r, c]) => grok.shell.tv.dataFrame.get(c, r - 1), [row, column] as
      [number, string]), {message: `row ${row} of ${column} changed`}).not.toBe(before);
  }, {description: 'the cell\'s value (row 1 the first) differs from the one named'});

export const valueContains = Then('the value in row {int} of {string} column contains {string}',
  async (page: Page, row: number, column: string, part: string) => {
    await expect.poll(() => page.evaluate(([r, c]) => String(grok.shell.tv.dataFrame.get(c, r - 1)),
      [row, column] as [number, string]), {message: `row ${row} of ${column}`}).toContain(part);
  }, {description: 'the cell\'s value (row 1 the first) holds the text named'});

/** In the page: whether Chem reads the cell as a molecule (Chem:isMolecule-like: the engine's get_mol, not a query). */
async function readsAsMolecule(page: Page, row: number, column: string): Promise<boolean> {
  return page.evaluate(async ([r, c]) => {
    const v = grok.shell.tv.dataFrame.get(c, r - 1);
    // canonical SMILES through Chem's canonicalizer (Chem:convertMolNotation refuses a SMILES to SMILES):
    // an unreadable or query-only value gives none
    const out = await grok.functions.call('Chem:canonicalize', {molecule: v}).catch(() => '');
    // query syntax in the SMILES itself (a CX extension after it holds commas of its own)
    return typeof out === 'string' && out !== '' && !/MALFORMED|malformed/.test(out) &&
      !/[,;!&~]|\[#|\$\(/.test(out.split(' |')[0]);
  }, [row, column] as [number, string]);
}

export const readsMolecule = Then('Chem reads row {int} of {string} as a molecule, not a query',
  async (page: Page, row: number, column: string) => {
    expect(await readsAsMolecule(page, row, column), `row ${row} of ${column} read as a molecule`).toBe(true);
  }, {description: 'Chem converts the cell to canonical SMILES (a query or an unreadable value gives none)'});

export const readsMolecule2 = Then('the value in row {int} of {string} column is read by Chem as a molecule',
  async (page: Page, row: number, column: string) => {
    expect(await readsAsMolecule(page, row, column), `row ${row} of ${column} read as a molecule`).toBe(true);
  }, {description: 'Chem converts the cell to canonical SMILES (a query or an unreadable value gives none)'});

// ---------------------------------------------------------------------------------------------------- Checks

/** The Checks section's lines for the current cell, once the context panel shows it: null when it has no Checks
 * section. The section is expanded when it is collapsed. */
/** In the page: the current cell made the context panel's object, as a click on it does; with no current cell (an
 * editor's OK can leave none: green2), the first chemical column's first row made current first. */
async function showCurrentCell(page: Page): Promise<void> {
  await page.evaluate(() => {
    const df = grok.shell.tv?.dataFrame;
    if (!df)
      return;
    let cell = df.currentCell;
    if (!cell?.column || cell.rowIndex < 0) {
      const col = df.columns.toList().find((c: any) => c.semType === DG.SEMTYPE.MOLECULE ||
        c.semType === DG.SEMTYPE.CHEMICAL_REACTION);
      if (!col)
        return;
      cell = df.cell(0, col.name);
      df.currentCell = cell;
    }
    grok.shell.o = DG.SemanticValue.fromTableCell(cell);
  });
}

async function checksLines(page: Page): Promise<{level: string, text: string}[] | null> {
  // the context panel shows the current cell. A click on a cell right after an editor's OK can leave the panel on the
  // value the editor had (seen after the sketcher dialog's OK): the cell is then made the panel's object, as the click
  // means to
  const shown = () => page.evaluate(() => {
    const df = grok.shell.tv?.dataFrame;
    const o = grok.shell.o;
    return !!df?.currentCell && df.currentCell.rowIndex >= 0 && o?.value === df.currentCell.value;
  });
  try {
    await expect.poll(shown, {message: 'the context panel shows the current cell', timeout: pollMs(5000)}).toBe(true);
  } catch {
    await showCurrentCell(page);
    await expect.poll(shown, {message: 'the context panel shows the current cell', timeout: pollMs(15000)}).toBe(true);
  }
  // the panel's sections are built a moment after: they settle (the same headers twice, 500 ms apart)
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')]
      .map((h) => h.textContent!.trim()).join('|'));
    const same = now === last && now !== '';
    last = now;
    return same;
  }, {message: 'the context panel\'s sections settled', intervals: [500], timeout: pollMs(20000)}).toBe(true);
  return page.evaluate(async () => {
    const header = [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')]
      .find((h) => h.textContent!.trim() === 'Checks') as HTMLElement | undefined;
    if (!header)
      return null;
    const pane = header.closest('.d4-accordion-pane')!;
    if (!pane.querySelector('.chem-footgun-checks'))
      header.click();
    for (let i = 0; i < 50 && !pane.querySelector('.chem-footgun-checks'); i++)
      await new Promise((r) => setTimeout(r, 100));
    return [...pane.querySelectorAll('.chem-footgun-checks > [data-level]')].map((l) => ({
      level: (l as HTMLElement).dataset.level!, text: l.textContent!.trim()}));
  });
}

const checksHave = async (page: Page, level: string, words: string[]) => {
  let lines: {level: string, text: string}[] | null = null;
  const message = `a ${level} mentioning ${words.map((w) => `"${w}"`).join(' and ')} in the Checks section`;
  await withState(expect.poll(async () => {
    lines = await checksLines(page);
    return !!lines?.some((l) => l.level === level && hasWords(l.text, words));
  }, {message, timeout: pollMs(20000)}).toBe(true), () => `${message}; it lists ${JSON.stringify(lines)}`);
};

export const checksWord = Then('the context panel\'s Checks section lists a {word} mentioning {string}',
  (page: Page, level: string, word: string) => checksHave(page, level, [word]),
{description: 'the Checks section of the current cell lists a flag of that level whose words hold the word (any case)'});

export const checksWords = Then('the context panel\'s Checks section lists a {word} mentioning {string} and {string}',
  (page: Page, level: string, a: string, b: string) => checksHave(page, level, [a, b]),
{description: 'the Checks section of the current cell lists a flag of that level whose words hold both words'});

export const checksForCell = Then('the context panel\'s Checks section for that cell lists a {word} mentioning {string}',
  async (page: Page, level: string, word: string) => {
    // the cell the editor just wrote: the context panel may still show the value it had, so the cell is made the
    // panel's object again, as a click on it does
    // a click on the cell, as the user makes to see its panel (the panel could keep the value the editor opened
    // with); then the cell made the panel's object, should the click land while the grid is still redrawing
    const column = await page.evaluate(() => {
      const df = grok.shell.tv?.dataFrame;
      const current = df?.currentCell;
      const name = current?.column && current.rowIndex >= 0 ? current.column.name : df?.columns.toList()
        .find((c) => c.semType === DG.SEMTYPE.MOLECULE || c.semType === DG.SEMTYPE.CHEMICAL_REACTION)?.name;
      return name ?? null;
    });
    if (column)
      await clickArea(page, `cell 1 of ${column}`, el('grid'));
    await showCurrentCell(page);
    await checksHave(page, level, [word]);
  },
{description: 'the cell the editor wrote, made the context panel\'s object again: its Checks section lists a flag of that level whose words hold the word'});

/** Crux mode: Chem's Checks exist (the section's two condition functions), so a cell with no Checks section is one the
 * checker found nothing to say about, not one no checker looked at (without this, the two steps below passed on a
 * build with no checker). */
async function checkerThere(page: Page): Promise<void> {
  if (await engineOf(page) !== 'Crux')
    return;
  const there = await page.evaluate(() => ['moleculeChecksShown', 'reactionChecksShown']
    .every((name) => DG.Func.find({package: 'Chem', name}).length > 0));
  expect(there, 'Chem\'s Checks in Crux mode (Chem:moleculeChecksShown and Chem:reactionChecksShown)').toBe(true);
}

export const checksNone = Then('the context panel\'s Checks section lists no flag', async (page: Page) => {
  await checkerThere(page);
  const lines = await checksLines(page);
  expect(lines ?? [], 'the Checks section lists nothing (or there is none)').toEqual([]);
}, {description: 'Crux mode: Chem\'s Checks exist; the current cell has no Checks section, or one that lists nothing'});

export const checksNoWord = Then('the context panel\'s Checks section lists no flag mentioning {string}',
  async (page: Page, word: string) => {
    await checkerThere(page);
    const lines = await checksLines(page);
    expect((lines ?? []).filter((l) => hasWords(l.text, [word])),
      `the Checks section lists no flag mentioning "${word}"; it lists ${JSON.stringify(lines)}`).toEqual([]);
  },
{description: 'Crux mode: Chem\'s Checks exist; the current cell\'s Checks section has no flag with the word'});

export const noChecks = Then('the context panel has no Checks section', async (page: Page) => {
  expect(await checksLines(page), 'no Checks section in the context panel').toBeNull();
}, {description: 'the context panel of the current cell, its sections settled, has no Checks section'});

// ---------------------------------------------------------------------------------------------------- reactions

/** The account's own reactions file (src/utils/reactions/reaction-storage.ts). */
const USER_REACTIONS = `(() => {
  const login = String(grok.shell.user.login).replace(/[^a-zA-Z0-9_-]/g, '_');
  return 'System:AppData/Chem/reactions/' + login + '-reactions.json';
})()`;

export const putBackReactions = Given('the user\'s saved reactions are put back when the feature ends', async (page: Page) => {
  const saved = await page.evaluate(async (file) => {
    const path = new Function(`return ${file}`)();
    return await grok.dapi.files.exists(path) ? await grok.dapi.files.readAsText(path) : null;
  }, USER_REACTIONS);
  atFeatureEnd(page, async () => {
    await page.evaluate(async ([file, text]) => {
      const path = new Function(`return ${file}`)();
      if (text === null) {
        if (await grok.dapi.files.exists(path))
          await grok.dapi.files.delete(path);
      } else
        await grok.dapi.files.writeAsText(path, text);
    }, [USER_REACTIONS, saved] as [string, string | null]);
  });
}, {tier: 'api', description: 'the account\'s saved reactions file as it is now, written back (or removed) when the feature ends'});

export const reactionSaved = Then('the reaction {string} is in the user\'s saved reactions', async (page: Page, name: string) => {
  // the account's file, or the open dialog's reaction browser, which reads the file again after a save: two workers
  // share the account, and the other's put-back can rewrite the file a moment later
  await expect.poll(() => page.evaluate(async ([file, n]) => {
    const card = [...document.querySelectorAll('.d4-dialog .d4-flex-col')].some((c) =>
      c.querySelector(':scope > canvas') && [...c.querySelectorAll(':scope > label')].map((l) => l.textContent?.trim())
        .includes(n) && [...c.querySelectorAll(':scope > label')].some((l) => l.textContent?.trim() === 'User'));
    if (card)
      return true;
    const path = new Function(`return ${file}`)();
    if (!await grok.dapi.files.exists(path))
      return false;
    return (JSON.parse(await grok.dapi.files.readAsText(path)) as any[]).some((r) => r.name === n);
  }, [USER_REACTIONS, name] as [string, string]), {message: `"${name}" saved`, timeout: pollMs(15000)}).toBe(true);
}, {description: 'the account\'s saved reactions file holds a reaction of that name, or the open dialog\'s reaction browser shows it as the user\'s'});

/** In the page: the reaction card of the topmost dialog's browser with that name, clicked. */
async function pickReaction(page: Page, name: string): Promise<void> {
  const card = page.locator('.d4-dialog').last().locator('.d4-flex-col:has(> canvas + label)', {
    has: page.locator(':scope > label', {hasText: new RegExp(`^\\s*${name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}\\s*$`)}),
  }).first();
  await expect(card, `the "${name}" reaction card`).toBeVisible({timeout: pollMs(30000)});
  await card.click();
}

/** Chem's reaction browser keeps the reactions it read for the page's life (reaction-storage.ts's cache; a save in its
 * New Reaction editor clears it), and the bdd harness keeps one page per worker: a reaction written to the account's
 * file shows only after such a save (else: "the "Three reactants" reaction card" not found). So a throwaway
 * reaction is saved through the dialog `menu` opens (the file, and so the throwaway, is put back when the feature ends),
 * and the dialog is closed. */
async function reloadSavedReactions(page: Page, menu: string): Promise<void> {
  await pickTopMenu(page, menu);
  const browser = page.locator('.d4-dialog').last();
  const add = browser.locator('button', {hasText: 'New Reaction'}).first();
  await expect(add, 'the dialog\'s New Reaction').toBeVisible({timeout: pollMs(30000)});
  await add.click();
  await fillEditorInput(page, 'Reaction SMARTS', '[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]');
  await fillEditorInput(page, 'Reaction Name', 'GUARD refresh');
  const editor = page.locator('.d4-dialog', {has: page.locator('.d4-dialog-title, .u2-dialog-title-text',
    {hasText: 'New Reaction'})}).last();
  await editor.locator('.d4-dialog-footer button', {hasText: /^\s*SAVE\s*$/}).click();
  await expect.poll(() => hasDialog(page, null, 'New Reaction'), {message: 'the New Reaction dialog closed',
    timeout: pollMs(30000)}).toBe(false);
  await page.locator('.d4-dialog').last().locator('.d4-dialog-footer button', {hasText: /^\s*CANCEL\s*$/}).click();
}

export const selectReaction = When('user selects the reaction {string} in the dialog', (page: Page, name: string) =>
  pickReaction(page, name),
{tier: 'ui', description: 'a click on the reaction card of that name in the topmost dialog\'s reaction browser'});

export const savedTwoComponent = Given('the user has a saved two-component reaction {string} {string}',
  async (page: Page, name: string, smarts: string) => {
    await page.evaluate(async ([file, n, s]) => {
      const path = new Function(`return ${file}`)();
      const list = await grok.dapi.files.exists(path) ? JSON.parse(await grok.dapi.files.readAsText(path)) : [];
      const id = `bdd-${n.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`;
      const reaction = {id, name: n, reactionSmarts: s, mode: 'two-component', category: 'Uncategorized',
        isUserDefined: true};
      await grok.dapi.files.writeAsText(path, JSON.stringify([...list.filter((r: any) => r.id !== id), reaction]));
    }, [USER_REACTIONS, name, smarts] as [string, string, string]);
    await reloadSavedReactions(page, 'Chem | Transform | Reactions | Two-Component Reaction...');
  }, {tier: 'api', description: 'a two-component reaction of that name and SMARTS among the account\'s saved reactions'});

/** In the page: every balloon's text from now on (window.__guardBalloons). */
async function watchBalloons(page: Page): Promise<void> {
  await page.evaluate(() => {
    const seen: string[] = [];
    (window as any).__guardBalloons = seen;
    document.querySelectorAll('.d4-balloon').forEach((b) => seen.push(b.textContent ?? ''));
    new MutationObserver((records) => {
      for (const r of records) {
        r.addedNodes.forEach((n) => {
          if (n instanceof HTMLElement && (n.matches('.d4-balloon') || n.querySelector('.d4-balloon')))
            seen.push(n.textContent ?? '');
        });
      }
    }).observe(document.body, {childList: true, subtree: true});
  });
}

export const runTransformation = When('user runs the transformation {string} on the {string} column',
  async (page: Page, name: string, column: string) => {
    await watchBalloons(page);
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await pickTopMenu(page, 'Chem | Transform | Reactions | Transformation...');
    await pickReaction(page, name);
    await page.evaluate((c) => {
      const d = [...document.querySelectorAll('.d4-dialog')].pop()!;
      const sel = [...d.querySelectorAll('select')].find((s) => [...(s as HTMLSelectElement).options]
        .some((o) => o.text === c)) as HTMLSelectElement | undefined;
      if (sel && sel.value !== c) {
        sel.value = c;
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }
    }, column);
    await page.evaluate((n) => (window as any).__guardColumnsBefore = n, before);
    await page.locator('.d4-dialog').last().locator('.d4-dialog-footer button[name="button-OK"]').click();
    // the run ends with a products column, or with a balloon (a run that fails: TPL-08 before its fix); what it made
    // is the Then steps' to say. The saved reactions' refresh says 'Reaction "GUARD refresh" saved.' a moment late,
    // which is no end of the run (a wait that ended on it ended before the run's column came)
    await expect.poll(() => page.evaluate((b) => grok.shell.tv.dataFrame.columns.length > b ||
      ((window as any).__guardBalloons ?? []).some((t: string) => !/Reaction ".*" saved/.test(t)), before),
    {message: 'the run ended: a products column added or a balloon shown', timeout: pollMs(60000)}).toBe(true);
  }, {tier: 'ui', description: 'Chem | Transform | Reactions | Transformation... on the column named, the reaction picked by its card, OK; the run ends with a products column added or a balloon shown (the balloons shown meanwhile are kept)'});

// ---------------------------------------------------------------------------------------------------- the Reaction pane

/** In the page: the Reaction pane of the context panel (expanded), its picture canvases' reactions in order: the
 * whole first, then each step's under Steps. */
async function paneReactions(page: Page): Promise<string[]> {
  return page.evaluate(async () => {
    const header = [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')]
      .find((h) => h.textContent!.trim() === 'Reaction') as HTMLElement | undefined;
    if (!header)
      return [];
    const pane = header.closest('.d4-accordion-pane')!;
    if (!pane.querySelector('canvas[data-reaction]'))
      header.click();
    for (let i = 0; i < 50 && !pane.querySelector('canvas[data-reaction]'); i++)
      await new Promise((r) => setTimeout(r, 100));
    return [...pane.querySelectorAll('canvas[data-reaction]')].map((c) => (c as HTMLElement).dataset.reaction!);
  });
}

export const paneSteps = Then('the Reaction panel shows {int} steps', async (page: Page, n: number) => {
  await expect.poll(async () => (await paneReactions(page)).length - 1,
    {message: 'the Reaction pane\'s steps', timeout: pollMs(20000)}).toBe(n);
}, {description: 'the Reaction pane draws the reaction and, under Steps, that many step pictures'});

export const paneStep = Then('step {int} of the Reaction panel draws {string}', async (page: Page, i: number,
  text: string) => {
  await expect.poll(async () => (await paneReactions(page))[i] ?? null,
    {message: `step ${i} of the Reaction pane`, timeout: pollMs(20000)}).toBe(text);
}, {description: 'the Reaction pane\'s step picture (step 1 the first under Steps) draws exactly that text'});

export const openFullScreen = When('user opens the Reaction panel\'s full-screen view', async (page: Page) => {
  await page.locator('.grok-prop-panel .d4-accordion-pane:has(canvas[data-reaction]) .fa-expand-alt').first().click();
  await expect(page.locator('.d4-dialog canvas[data-reaction]').first(), 'the full-screen view').toBeVisible();
}, {tier: 'ui', description: 'a click on the Reaction pane\'s full-screen icon'});

export const fullScreenStep = Then('the full-screen view\'s step {int} draws {string}', async (page: Page, i: number,
  text: string) => {
  await expect.poll(() => page.evaluate((k) => {
    const d = [...document.querySelectorAll('.d4-dialog')].filter((x) => x.querySelector('canvas[data-reaction]')).pop();
    return d ? ([...d.querySelectorAll('canvas[data-reaction]')][k] as HTMLElement | undefined)?.dataset.reaction ??
      null : null;
  }, i), {message: `step ${i} of the full-screen view`}).toBe(text);
}, {description: 'the full-screen view\'s step picture (step 1 the first under Steps) draws exactly that text'});

export const noReactionPane = Then('the context panel has no Reaction panel', async (page: Page) => {
  await checksLines(page); // the panel shows the current cell, its sections settled
  const headers = await page.evaluate(() => [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')]
    .map((h) => h.textContent!.trim()));
  expect(headers, 'the context panel\'s sections').not.toContain('Reaction');
}, {description: 'the context panel of the current cell, its sections settled, has no Reaction section'});

// ---------------------------------------------------------------------------------------------------- the enumerator
/* The Reaction Enumerator (Chem | Reactions): its input tables are table pickers (`select[name="input-…"]`) listing
   the workspace's tables, so a step adds its table to the workspace (the view stays the enumerator's) and picks it, as
   a user picks an open table; the column pickers follow by the config's column names (reaction_smarts, blocking_fg,
   reaction_name; SMILES; SMARTS). The words beside Enumerate are `.chem-footgun-enumerate` in the ribbon (its whole
   text in `aria-label`, its level in `data-level`), as a probe of a running Datagrok found them. */

/** The picker of each input table, and the ribbon chip that opens its pane when it is not shown. */
const ENUMERATOR_INPUTS: {[what: string]: {select: string, chip: string}} = {
  templates: {select: 'input-Reaction-templates-file', chip: 'reactions'},
  buildingBlocks: {select: 'input-Building-blocks-file', chip: 'BBs'},
  exclusion: {select: 'input-Exclusion-SMARTS-file-(optional)', chip: 'Reagents'},
};

/** The enumerator's ribbon Enumerate button. */
const enumerateButton = (page: Page) => page.locator('.d4-ribbon-item button', {hasText: /^\s*Enumerate\s*$/}).first();

/** The words beside Enumerate: whether they are there, their level and their whole text. */
async function enumerateWords(page: Page): Promise<{present: boolean, level: string, whole: string}> {
  return page.evaluate(() => {
    const words = [...document.querySelectorAll('.chem-footgun-enumerate')]
      .find((e) => (e as HTMLElement).offsetParent !== null) as HTMLElement | undefined;
    return words ? {present: true, level: words.dataset.level ?? '', whole: words.getAttribute('aria-label') ?? ''} :
      {present: false, level: '', whole: ''};
  });
}

let enumeratorTables = 0;

/** Adds a table of `columns` to the workspace, the enumerator's view kept current, and picks it in the picker of
 * `what`; resolves once the picker holds it. */
async function pickEnumeratorTable(page: Page, what: string, columns: {[name: string]: string[]}): Promise<void> {
  const name = `guard-${what}-${++enumeratorTables}`;
  await page.evaluate(([n, cols]) => {
    const view = grok.shell.v;
    const df = DG.DataFrame.fromColumns(Object.entries(cols).map(([c, values]) => DG.Column.fromStrings(c, values)));
    df.name = n;
    grok.shell.addTable(df);
    if (grok.shell.v !== view && view)
      grok.shell.v = view;
  }, [name, columns] as [string, {[name: string]: string[]}]);
  const {select, chip} = ENUMERATOR_INPUTS[what];
  const picker = page.locator(`.chem-enumerator select[name="${select}"]`).first();
  if (!await picker.isVisible()) {
    await page.locator('.d4-ribbon-item', {hasText: new RegExp(chip)}).first().click();
    try {
      await expect(picker).toBeVisible({timeout: pollMs(5000)});
    } catch {
      // the chip's pane is open, and the picker sits in a collapsed accordion pane inside it: open those, outermost first
      await picker.evaluate((el) => {
        const panes: HTMLElement[] = [];
        for (let p = el.closest('.d4-accordion-pane') as HTMLElement | null; p;
          p = p.parentElement?.closest('.d4-accordion-pane') as HTMLElement | null)
          panes.unshift(p);
        for (const p of panes) {
          const content = p.querySelector(':scope > .d4-accordion-pane-content') as HTMLElement | null;
          const header = p.querySelector(':scope > .d4-accordion-pane-header') as HTMLElement | null;
          if (header && (!content || content.getClientRects().length === 0))
            header.click();
        }
      });
    }
    await expect(picker, `the ${what} picker`).toBeVisible({timeout: pollMs(10000)});
  }
  await expect.poll(async () => (await picker.locator('option').allTextContents()).includes(name),
    {message: `"${name}" among the ${what} picker's tables`, timeout: pollMs(10000)}).toBe(true);
  await picker.selectOption({label: name});
  await expect(picker, `the ${what} picker holds "${name}"`).toHaveValue(/.+/, {timeout: pollMs(5000)});
  await expect.poll(() => picker.evaluate((s: HTMLSelectElement) => s.options[s.selectedIndex]?.label),
    {message: `the ${what} picker holds "${name}"`, timeout: pollMs(10000)}).toBe(name);
}

export const openEnumerator = Given('user opens the Reaction Enumerator from {string}', async (page: Page, path: string) => {
  const browsePath = await page.evaluate(async () => {
    const f = DG.Func.find({package: 'Chem', name: 'reactionEnumeratorApp'})[0];
    const view = await f.apply();
    if (view)
      grok.shell.addView(view);
    return f.options?.browsePath ?? '';
  });
  expect(browsePath, 'the Reaction Enumerator\'s place in Browse').toBe(path);
  // the app opens on "How to combine": the templates picker is there, in a pane its ribbon chip opens
  await expect(page.locator(`.chem-enumerator select[name="${ENUMERATOR_INPUTS.templates.select}"]`).first(),
    'the Reaction Enumerator, its templates picker').toBeAttached({timeout: pollMs(60000)});
  await expect(enumerateButton(page), 'the Enumerate button').toBeVisible({timeout: pollMs(30000)});
}, {tier: 'api', description: 'the Reaction Enumerator app (its Browse path checked) opened as its app function opens it, with its default files, its templates picker and Enumerate shown'});

export const enumeratorTemplates = Given('the enumerator\'s templates are:', async (page: Page, rows: string[][]) => {
  const [header, ...body] = rows;
  const columns: {[name: string]: string[]} = {};
  header.forEach((h, i) => columns[h.trim()] = body.map((r) => (r[i] ?? '').trim()));
  await pickEnumeratorTable(page, 'templates', columns);
}, {tier: 'ui', description: 'a table of the rows (the header names its columns: reaction_smarts, blocking_fg, reaction_name, the config\'s names) added to the workspace and picked as the templates file'});

export const enumeratorBuildingBlocks = Given('the enumerator\'s building blocks are {string}',
  async (page: Page, list: string) => {
    await pickEnumeratorTable(page, 'buildingBlocks', {SMILES: list.split(',').map((s) => s.trim()).filter((s) => s)});
  }, {tier: 'ui', description: 'a table of one SMILES column holding the comma-separated building blocks, added to the workspace and picked as the building blocks file'});

export const enumeratorExclusion = Given('the enumerator\'s exclusion SMARTS are {string}',
  async (page: Page, list: string) => {
    await pickEnumeratorTable(page, 'exclusion', {SMARTS: list.split(', ').map((s) => s.trim()).filter((s) => s)});
  }, {tier: 'ui', description: 'a table of one SMARTS column holding the SMARTS (separated by ", ": a SMARTS may hold a comma), added to the workspace and picked as the exclusion SMARTS file'});

export const enumerateDisabled = Then('the Enumerate button is disabled', async (page: Page) => {
  await expect.poll(() => enumerateButton(page).evaluate((b: HTMLButtonElement) => b.disabled),
    {message: 'Enumerate disabled', timeout: pollMs(15000)}).toBe(true);
}, {description: 'the ribbon\'s Enumerate button is disabled (a click runs nothing)'});

export const enumerateEnabled = Then('the Enumerate button is enabled', async (page: Page) => {
  await expect.poll(() => enumerateButton(page).evaluate((b: HTMLButtonElement) => b.disabled),
    {message: 'Enumerate enabled', timeout: pollMs(15000)}).toBe(false);
}, {description: 'the ribbon\'s Enumerate button is enabled'});

export const enumerateWordsUnreadable = Then('the words beside the Enumerate button name {string} and say it cannot be read',
  async (page: Page, smarts: string) => {
    await expect.poll(async () => {
      const w = await enumerateWords(page);
      return w.present && w.level === 'block' && w.whole.includes(smarts) && /cannot be read/i.test(w.whole);
    }, {message: `the words beside Enumerate name ${smarts} as unreadable (a block)`, timeout: pollMs(15000)}).toBe(true);
  }, {description: 'the words beside Enumerate are a block whose whole text names the SMARTS, as written, and says it cannot be read'});

export const enumerateWordsNotName = Then('the words beside the Enumerate button do not name {string}',
  async (page: Page, smarts: string) => {
    const w = await enumerateWords(page);
    expect(w.present, 'the words beside Enumerate (Crux mode)').toBe(true);
    expect(w.whole, 'the words beside Enumerate').not.toContain(smarts);
  }, {description: 'the words beside Enumerate are there and their whole text does not hold the SMARTS'});

export const enumerateWordsWarnNone = Then(
  'the words beside the Enumerate button warn that {string} matches none of the building blocks',
  async (page: Page, smarts: string) => {
    await expect.poll(async () => {
      const w = await enumerateWords(page);
      return w.present && w.level === 'warn' && w.whole.includes(smarts) &&
        /matches none of the building blocks/i.test(w.whole);
    }, {message: `the words beside Enumerate warn that ${smarts} matches no building block`, timeout: pollMs(15000)})
      .toBe(true);
  }, {description: 'the words beside Enumerate are a warn whose whole text names the SMARTS and says it matches none of the building blocks'});

export const enumerateNoWords = Then('no foot-gun words show beside the Enumerate button', async (page: Page) => {
  await expect(enumerateButton(page), 'the Enumerate button').toBeVisible({timeout: pollMs(15000)});
  const w = await enumerateWords(page);
  expect(w.present && w.whole !== '' ? w.whole : '', 'the words beside Enumerate').toBe('');
}, {description: 'beside Enumerate there are no foot-gun words (none, or empty)'});

export const changeBlocking = When('user changes the blocking SMARTS of template {int} to {string}',
  async (page: Page, row: number, smarts: string) => {
    await page.evaluate(([r, v, select]) => {
      const picker = document.querySelector(`.chem-enumerator select[name="${select}"]`) as HTMLSelectElement;
      const name = picker.options[picker.selectedIndex].label;
      const df = (grok.shell.tables as any[]).find((t) => t.name === name);
      df.set('blocking_fg', r - 1, v);
    }, [row, smarts, ENUMERATOR_INPUTS.templates.select] as [number, string, string]);
  }, {tier: 'api', description: 'the blocking_fg cell of that template row (1-based) of the picked templates table set to the SMARTS, as an edit of the cell sets it'});

// ---------------------------------------------------------------------------------------------------- the reaction cell editor

/* The reaction cell editor: a double-click on a reaction cell opens, in either mode, Crux Sketch's reaction
   document in a dialog. Its box `.chem-reaction-cell-editor` is a widget named cruxSketch, whose status gives an area
   "atom N" per atom (by index, in the widget's own px) and the readings ready, pending, atoms, bonds, arrows and changes
   (a person's edits). Crux Sketch's controls go by their test ids in its shadow root: the eraser `toolbar.erase`, Undo
   `toolbar.undo`, the single bond tool `toolbar.bond.single`; an arrow is a `g[data-arrow]` group (a click on its
   box's centre lands on a straight arrow's shaft). Why Crux Sketch did not open a cell's reaction is the
   line `.chem-reaction-cell-editor-refused` above its canvas. */

/** In the page: the reaction editor's readings (its widget's status values), or null while none is open. */
async function reactionEditorReadings(page: Page): Promise<{[name: string]: any} | null> {
  return page.evaluate((sel) => {
    const dialogs = [...document.querySelectorAll(sel)];
    const root = dialogs[dialogs.length - 1]?.querySelector('.chem-reaction-cell-editor');
    return root ? DG.Widget.find(root)?.getWidgetStatus()?.values ?? null : null;
  }, REACTION_EDITOR);
}

/** Waits for the reaction editor to have read the cell and to hold still (ready, not pending); its readings. */
async function reactionEditorReady(page: Page): Promise<{[name: string]: any}> {
  let r: {[name: string]: any} | null = null;
  await withState(expect.poll(async () => {
    r = await reactionEditorReadings(page);
    return r !== null && r.ready === true && r.pending === false;
  }, {message: 'the reaction editor ready', timeout: pollMs(30000)}).toBe(true), () => `readings: ${JSON.stringify(r)}`);
  return r!;
}

/** In Crux mode, waits for the reaction editor to be ready; in RDKit mode (a text dialog), for its dialog. */
async function reactionEditorOpen(page: Page): Promise<void> {
  await expect(page.locator(REACTION_EDITOR).last(), 'the reaction editor').toBeVisible({timeout: pollMs(30000)});
  if (await engineOf(page) === 'Crux')
    await reactionEditorReady(page);
}

/** In the page: the centres of the reaction editor's atoms (by index; null for one not drawn) and its arrows' boxes, in
 * the page's CSS px. */
async function editorGeometry(page: Page): Promise<{atoms: ({x: number, y: number} | null)[], arrows: Rect[]}> {
  return page.evaluate((sel) => {
    const dialogs = [...document.querySelectorAll(sel)];
    const root = dialogs[dialogs.length - 1]?.querySelector('.chem-reaction-cell-editor') as HTMLElement | undefined;
    if (!root)
      throw new Error('no reaction editor is open');
    const status = DG.Widget.find(root).getWidgetStatus();
    const at = root.getBoundingClientRect();
    const atoms: ({x: number, y: number} | null)[] = [];
    for (let i = 0; i < (status.values?.atoms ?? 0); i++) {
      const a = status.hitAreas[`atom ${i}`];
      atoms.push(a ? {x: at.left + a.x + a.width / 2, y: at.top + a.y + a.height / 2} : null);
    }
    const sketch = root.querySelector('crux-sketch');
    const arrows = [...(sketch?.shadowRoot?.querySelectorAll('g[data-arrow]') ?? [])].map((g) => {
      const r = g.getBoundingClientRect();
      return {left: r.left, top: r.top, width: r.width, height: r.height};
    });
    return {atoms, arrows};
  }, REACTION_EDITOR);
}

/** A click on one of Crux Sketch's controls in the reaction editor, by its test id. */
const editorControl = (page: Page, id: string) =>
  page.locator(`${REACTION_EDITOR} crux-sketch [data-testid="${id}"]`).last().click();

/** Waits for a person's edit to have reached the reaction editor (its `changes` past `before`), then for it to settle. */
async function editorChanged(page: Page, before: number, what: string): Promise<void> {
  await expect.poll(async () => (await reactionEditorReadings(page))?.changes ?? -1,
    {message: `${what}: a change of the person's`, timeout: pollMs(10000)}).toBeGreaterThan(before);
  await reactionEditorReady(page);
}

/** The atoms of reaction text (SMILES or SMARTS; a CX extension left out) in the order they are written, which is the
 * order the reaction editor's atoms come in (checked by side wherever a step uses it): each its role, component,
 * element and map number. A bracket atom's element is its first symbol (`[NX3;H2,H1;...]` N, `[#6:2]` C); recursive
 * SMARTS inside it are skipped whole. */
function writtenAtoms(text: string): {role: string, component: number, element: string, map: number | null}[] {
  const parts = text.split(' |')[0].trim().split('>');
  if (parts.length !== 3)
    throw new Error(`not one reaction step: ${text}`);
  const numbers: {[n: string]: string} = {5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br'};
  const out: {role: string, component: number, element: string, map: number | null}[] = [];
  parts.forEach((part, p) => {
    const role = ['reactant', 'agent', 'product'][p];
    let component = 0;
    for (let i = 0; i < part.length;) {
      if (part[i] === '[') {
        let depth = 1;
        let j = i + 1;
        for (; j < part.length && depth > 0; j++)
          depth += part[j] === '[' ? 1 : part[j] === ']' ? -1 : 0;
        const inner = part.slice(i + 1, j - 1);
        const symbol = /^\d*(#(\d+)|Cl|Br|[A-Z][a-z]?|[a-z])/.exec(inner);
        const element = !symbol ? '?' : symbol[2] ? numbers[symbol[2]] ?? `#${symbol[2]}` :
          symbol[1][0].toUpperCase() + symbol[1].slice(1);
        const map = /:(\d+)$/.exec(inner);
        out.push({role, component, element, map: map ? +map[1] : null});
        i = j;
        continue;
      }
      if (part[i] === '.') {
        component++;
        i++;
        continue;
      }
      const organic = /^(Cl|Br|[BCNOPSFI]|[bcnops]|\*)/.exec(part.slice(i));
      if (organic) {
        out.push({role, component, element: organic[1].toUpperCase(), map: null});
        i += organic[1].length;
        continue;
      }
      i++;
    }
  });
  return out;
}

/** The value of row 1 of the table the scenario opened (the cell the editor was opened on). */
function openedValue(page: Page): string {
  const columns = opened.get(page);
  const first = columns ? Object.values(columns)[0]?.[0] : undefined;
  if (first === undefined)
    throw new Error('no table was opened earlier in the scenario');
  return first;
}

/** Adds a methyl to atom `index` of the reaction editor as a person does: the single bond tool, a click on the atom
 * (crux-sketch BOND-003: a bond to a new carbon). `side` is where the atom must be: left of the arrow (a reactant's) or
 * right of it (a product's). */
async function addMethyl(page: Page, index: number, side: 'reactant' | 'product'): Promise<void> {
  const before = await reactionEditorReady(page);
  const {atoms, arrows} = await editorGeometry(page);
  const a = atoms[index];
  expect(a, `atom ${index} drawn (atoms: ${atoms.length})`).toBeTruthy();
  expect(arrows.length, 'one arrow drawn').toBe(1);
  const arrow = arrows[0];
  if (side === 'reactant')
    expect(a!.x, `atom ${index} left of the arrow (a reactant's)`).toBeLessThan(arrow.left);
  else
    expect(a!.x, `atom ${index} right of the arrow (a product's)`).toBeGreaterThan(arrow.left + arrow.width);
  await editorControl(page, 'toolbar.bond.single');
  await page.mouse.move(a!.x - 3, a!.y - 3);
  await page.mouse.move(a!.x, a!.y);
  await page.mouse.click(a!.x, a!.y);
  await editorChanged(page, before.changes, `a methyl on atom ${index}`);
  expect((await reactionEditorReadings(page))?.atoms, 'one atom more').toBe(before.atoms + 1);
}

export const reactionEditorSlotEmpty = Then('the footer of the reaction editor has an empty middle slot',
  async (page: Page) => {
    await reactionEditorOpen(page);
    await expectSlot(page, 'reaction', undefined, 'an empty middle slot in the reaction editor\'s footer',
      (s) => s.present && s.text === '');
  }, {description: 'once the reaction editor is ready, its footer holds the middle slot, and it says nothing'});

export const noteReactionEditorButtons = When('user notes where the {word} and {word} buttons of the reaction editor are',
  async (page: Page, a: string, b: string) => {
    await reactionEditorOpen(page);
    await noteButtons(page, a, b, 'reaction');
  }, {tier: 'api', description: 'the boxes of the two footer buttons of the reaction editor, once it is ready and they hold still, kept for a later step'});

export const reactionEditorButtonsStill = Then('the {word} and {word} buttons of the reaction editor are where they were',
  (page: Page, a: string, b: string) => buttonsStill(page, a, b, 'reaction'),
{description: 'each button\'s box (left, top, width, height) within 0.5 px of the one noted'});

const reactionEditorButtonIs = (want: boolean) => async (page: Page, name: string) => {
  await reactionEditorOpen(page);
  await expect.poll(() => buttonDisabled(page, name, 'reaction'),
    {message: `${name} of the reaction editor ${want ? 'disabled' : 'enabled'}`}).toBe(want);
};
export const reactionEditorButtonDisabled = Then('the {word} button of the reaction editor is disabled',
  reactionEditorButtonIs(true), {description: 'the reaction editor\'s footer button is disabled (a disabled OK neither clicks nor answers Enter)'});
export const reactionEditorButtonEnabled = Then('the {word} button of the reaction editor is enabled',
  reactionEditorButtonIs(false), {description: 'the reaction editor\'s footer button is enabled'});

export const clickReactionEditorButton = When('user clicks on {word} button of the reaction editor',
  async (page: Page, name: string) => {
    await reactionEditorOpen(page);
    await page.locator(REACTION_EDITOR).last().locator('.d4-dialog-footer button',
      {hasText: new RegExp(`^\\s*${name}\\s*$`, 'i')}).click();
    await expect(page.locator(REACTION_EDITOR), 'the reaction editor closed').toHaveCount(0, {timeout: pollMs(10000)});
  }, {tier: 'ui', description: 'once the reaction editor is ready, a click on its footer button (OK, CANCEL); the dialog closes'});

export const eraseArrow = When('user erases the arrow in the reaction editor', async (page: Page) => {
  const before = await reactionEditorReady(page);
  const {arrows} = await editorGeometry(page);
  expect(arrows.length, 'an arrow drawn in the reaction editor').toBeGreaterThan(0);
  await editorControl(page, 'toolbar.erase');
  const a = arrows[0];
  await page.mouse.click(a.left + a.width / 2, a.top + a.height / 2);
  await editorChanged(page, before.changes, 'the arrow erased');
  expect((await reactionEditorReadings(page))?.arrows, 'arrows drawn').toBe(arrows.length - 1);
}, {tier: 'ui', description: 'Crux Sketch\'s eraser, then a click on the middle of the arrow\'s box (its shaft): the arrow erased, one edit'});

export const undoEdit = When('user undoes the last edit in the reaction editor', async (page: Page) => {
  const before = await reactionEditorReady(page);
  await editorControl(page, 'toolbar.undo');
  await editorChanged(page, before.changes, 'the last edit undone');
}, {tier: 'ui', description: 'a click on Crux Sketch\'s Undo: the last edit undone, as a change of the person\'s'});

export const methylOnProductMap = When('user adds a methyl to the product\'s carbon {int} in the reaction editor',
  async (page: Page, map: number) => {
    const atoms = writtenAtoms(openedValue(page));
    const index = atoms.findIndex((a) => a.role === 'product' && a.map === map);
    expect(index, `a product atom with map number ${map} in ${openedValue(page)}`).toBeGreaterThanOrEqual(0);
    await addMethyl(page, index, 'product');
  }, {tier: 'ui', description: 'the single bond tool, then a click on the product\'s atom with that map number (found in the cell\'s text, its place checked right of the arrow): a methyl on it'});

export const methylOnFirstReactant = When('user adds a methyl to the first reactant in the reaction editor',
  (page: Page) => addMethyl(page, 0, 'reactant'),
{tier: 'ui', description: 'the single bond tool, then a click on atom 0 (its place checked left of the arrow): a methyl on the first reactant'});

export const noReactantFlag = Then('the footer shows no flag about the reactants', async (page: Page) => {
  await reactionEditorReady(page);
  await expectSlot(page, 'reaction', undefined, 'no flag about the reactants (no query feature lost)',
    (s) => s.present && !linesOf(s).some((l) => /query feature|reactant/i.test(l)));
}, {description: 'the reaction editor\'s slot is there, and no flag of it names a reactant or a query feature lost'});

export const valueStillContains = Then('the value in row {int} of {string} column still contains {string}',
  async (page: Page, row: number, column: string, part: string) => {
    await expect.poll(() => cellValue(page, row, column), {message: `row ${row} of ${column}`}).toContain(part);
  }, {description: 'the cell\'s value (row 1 the first) still holds the text named, as written'});

// a route Crux Sketch does not open (a step with no product; a convergent route opens since crux-core 15b1bbd65)
export const editorSaysRefused = Then('the reaction editor says it cannot open the route', async (page: Page) => {
  const line = page.locator(`${REACTION_EDITOR} .chem-reaction-cell-editor-refused`).last();
  await expect(line, 'the editor\'s one message, in place of its sketcher').toBeVisible({timeout: pollMs(30000)});
  const text = (await line.innerText()).toLowerCase();
  expect(/does not open|cannot (edit|open)/.test(text), `words saying it does not open it: ${text}`).toBe(true);
  expect(text, 'words naming the route').toContain('route');
}, {description: 'the editor\'s one message shows in place of its sketcher, saying it does not open the reaction and ' +
  'naming the route'});

/** The rows of CSV text (quoted fields with commas, doubled quotes and line breaks; a line break in a field read as
 * "\n"). */
function csvRows(text: string): string[][] {
  const rows: string[][] = [];
  let row: string[] = [];
  let field = '';
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (ch === '"')
        quoted = false;
      else if (ch !== '\r')
        field += ch;
    } else if (ch === '"')
      quoted = true;
    else if (ch === ',') {
      row.push(field);
      field = '';
    } else if (ch === '\n') {
      rows.push([...row, field]);
      row = [];
      field = '';
    } else if (ch !== '\r')
      field += ch;
  }
  if (field !== '' || row.length)
    rows.push([...row, field]);
  return rows;
}

/** Chem's reaction corpus (its package tests' file): id, kind, reaction. */
const REACTION_CORPUS = '../files/tests/reactions-crux.csv';

export const corpusRxnTable = Given('user opens a table {string} with a reaction column {string} of units {string} ' +
  'holding the corpus\'s {string} RXN block',
async (page: Page, name: string, column: string, units: string, id: string) => {
  const {readFileSync} = await import('node:fs');
  const [header, ...rows] = csvRows(readFileSync(REACTION_CORPUS, 'utf8'));
  const row = rows.find((r) => r[header.indexOf('id')] === id);
  if (!row)
    throw new Error(`no row "${id}" in ${REACTION_CORPUS}`);
  const value = row[header.indexOf('reaction')];
  expect(value.startsWith('$RXN'), `the corpus's "${id}" is an RXN block`).toBe(true);
  await openReactionTable(page, name, column, units, [value]);
}, {tier: 'api', description: 'a one-row reaction column of the units named holding the RXN block of Chem\'s reaction corpus (files/tests/reactions-crux.csv) by its id'});

export const rxnNotationWarn = Then('the middle slot of the footer shows a warn flag mentioning {string} and the ' +
  'notation it writes instead', async (page: Page, word: string) => {
  await expectSlot(page, 'reaction', undefined, `a warn mentioning "${word}" and the notation OK writes instead`,
    (s) => s.present && s.level === 'warn' &&
      linesOf(s).some((l) => hasWords(l, [word]) && /reaction (SMILES|SMARTS)/i.test(l)));
},
{description: 'the reaction editor\'s slot shows a warn (its worst) one of whose flags holds the word and names the notation OK writes instead (reaction SMILES or SMARTS)'});

// ---------------------------------------------------------------------------------------------------- the sketcher

export const slotShowsWordsNoFlag = Then('the middle slot of the footer shows a {word} mentioning {string} and {string}',
  (page: Page, level: string, a: string, b: string) => expectSlot(page, null, undefined,
    `a ${level} mentioning "${a}" and "${b}"`, (s) => s.present && linesOf(s).length > 0 &&
      s.level !== '' && linesOf(s).some((l) => hasWords(l, [a, b])) && slotHasLevel(s, level)),
{description: 'the middle slot holds a flag of that level whose words hold both words (any case)'});

/** Whether the slot holds a flag of `level` (its worst is shown; the rest are in its whole text, in level order). */
function slotHasLevel(s: SlotState, level: string): boolean {
  const order = ['block', 'warn', 'note'];
  return order.indexOf(s.level) <= order.indexOf(level);
}

export const conformerTable = Given(
  'user opens a table {string} with a molblock column {string} holding a 3D molfile of {string}',
  async (page: Page, name: string, column: string, smiles: string) => {
    await page.evaluate(([n, c, s]) => {
      // the molecule's 2D molblock, its z coordinates set apart from 0 (a conformer's)
      const flat: string = grok.chem.convert(s, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
      const lines = flat.split('\n');
      const atoms = +lines[3].slice(0, 3);
      for (let i = 4; i < 4 + atoms; i++)
        lines[i] = lines[i].slice(0, 20) + (0.5 * (i - 3)).toFixed(4).padStart(10) + lines[i].slice(30);
      const col = DG.Column.fromStrings(c, [lines.join('\n')]);
      col.semType = DG.SEMTYPE.MOLECULE;
      col.meta.units = DG.chem.Notation.MolBlock;
      const df = DG.DataFrame.fromColumns([col]);
      df.name = n;
      grok.shell.addTableView(df);
    }, [name, column, smiles] as [string, string, string]);
    await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid, name),
      {message: `a table view of "${name}"`}).toBe(true);
  }, {tier: 'api', description: 'a one-row molblock column holding the molecule\'s V2000 molblock with non-zero z coordinates (0.5 Å apart, atom by atom)'});

/** In the page: the cell's value (row 1 the first). */
const cellValue = (page: Page, row: number, column: string) =>
  page.evaluate(([r, c]) => String(grok.shell.tv.dataFrame.get(c, r - 1) ?? ''), [row, column] as [number, string]);

/** Whether a molblock's atoms have a z coordinate other than 0. */
function hasZ(molblock: string): boolean {
  if (molblock.includes('V3000'))
    return [...molblock.matchAll(/^M {2}V30 \d+ \S+ \S+ \S+ (\S+)/gm)].some((m) => Math.abs(+m[1]) > 1e-4);
  const lines = molblock.split('\n');
  const atoms = +(lines[3] ?? '').slice(0, 3);
  return lines.slice(4, 4 + atoms).some((l) => Math.abs(+l.slice(20, 30)) > 1e-4);
}

export const keeps3D = Then('the cell in row {int} of {string} still holds its 3D coordinates',
  async (page: Page, row: number, column: string) => {
    expect(hasZ(await cellValue(page, row, column)), 'a z coordinate other than 0 in the molblock').toBe(true);
  }, {description: 'the cell\'s molblock has an atom whose z is not 0'});

export const atomList = Then('the cell in row {int} of {string} holds a molblock with an atom list',
  async (page: Page, row: number, column: string) => {
    const v = await cellValue(page, row, column);
    expect(v.includes('M  END') && (/^M {2}ALS/m.test(v) || /\[[A-Z][a-z]?,/.test(v) || / L {3}0/.test(v)),
      `an atom list in ${JSON.stringify(v.slice(0, 200))}...`).toBe(true);
  }, {description: 'the cell holds a molblock with an atom list (V2000 M  ALS, or a V3000 list atom)'});

export const putChain = When('user puts a chain of {int},{int} carbon atoms into the Crux sketcher of the dialog',
  (page: Page, thousands: number, rest: number) => putIntoCrux(page, 'C'.repeat(thousands * 1000 + rest), false),
{tier: 'api', description: 'a linear alkane of that many carbons set on the sketcher dialog\'s host as SMILES'});

export const v3000Of = Then('the cell in row {int} of {string} holds a V3000 molblock of {int},{int} atoms',
  async (page: Page, row: number, column: string, thousands: number, rest: number) => {
    const v = await cellValue(page, row, column);
    const n = thousands * 1000 + rest;
    expect(v.includes('V3000'), 'a V3000 molblock').toBe(true);
    expect(+(/M {2}V30 COUNTS (\d+)/.exec(v)?.[1] ?? -1), 'its atom count').toBe(n);
  }, {description: 'the cell holds a V3000 molblock whose COUNTS line names that many atoms'});

// ---------------------------------------------------------------------------------------------------- mapping

/** In the page: the Reaction pane's Mapping section: its status lines and its flag lines (`[data-level]`). */
async function mappingState(page: Page): Promise<{status: string, flags: {level: string, text: string}[],
  state: string | null}> {
  return page.evaluate(() => {
    const header = [...document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')]
      .find((h) => h.textContent!.trim() === 'Reaction') as HTMLElement | undefined;
    const pane = header?.closest('.d4-accordion-pane');
    if (!pane)
      return {status: '', flags: [], state: null};
    const h3 = [...pane.querySelectorAll('h3')].find((h) => h.textContent!.trim() === 'Mapping');
    const section: Element[] = [];
    for (let e = h3?.nextElementSibling ?? null; e !== null && e.tagName !== 'H3'; e = e.nextElementSibling)
      section.push(e);
    // the automatic mapping: its status's state, waiting, mapping, done, timeout, failed or refused; null on a panel
    // with the Map atoms button
    const status = pane.querySelector('.chem-reaction-mapping') as HTMLElement | null;
    return {status: section.map((e) => e.textContent!.trim()).join('\n'),
      flags: section.flatMap((e) => [...e.querySelectorAll('[data-level]')]).map((l) => ({
        level: (l as HTMLElement).dataset.level!, text: l.textContent!.trim()})),
      state: status?.dataset.state ?? null};
  });
}

/** Expands the Reaction pane (once the panel shows the current cell). */
async function reactionPane(page: Page): Promise<void> {
  await checksLines(page);
  await paneReactions(page);
}

export const mapAtoms = When('user clicks the {string} button of the Reaction panel', async (page: Page, text: string) => {
  await reactionPane(page);
  await page.locator('.grok-prop-panel .d4-accordion-pane:has(canvas[data-reaction]) button', {hasText: text}).first()
    .click();
}, {tier: 'ui', description: 'a click on the Reaction pane\'s button of that text (the pane expanded first)'});

/** Whether the mapping has answered: the automatic mapping's state past waiting and mapping, or (a panel with the Map
 * atoms button) its words. */
const mappedDone = (s: {status: string, state: string | null}) => s.state !== null ?
  // the answers only: a state before one (waiting, loading the mapper, mapping) is no answer
  ['done', 'timeout', 'failed', 'refused'].includes(s.state) :
  /\bmapped\b|not mapped|Not mapped/i.test(s.status) && !/Mapping\.\.\./.test(s.status);

export const mappingAnswers = When('the Reaction panel\'s mapping answers', async (page: Page) => {
  await reactionPane(page);
  let s: {status: string, state: string | null} | null = null;
  await withState(expect.poll(async () => {
    s = await mappingState(page);
    return s.state !== null && mappedDone(s);
  }, {message: 'the Reaction panel\'s automatic mapping answered', timeout: pollMs(180000)}).toBe(true),
  () => `the Mapping section: ${JSON.stringify(s)}`);
}, {tier: 'api', description: 'the Reaction pane expanded, its automatic mapping (Crux mode) has answered: its status one of the answers (done, timeout, failed, refused), never waiting, loading the mapper or mapping'});

export const mappingConfidence = Then('the Reaction panel\'s Mapping section shows the confidence {string}',
  async (page: Page, confidence: string) => {
    // the status says how sure the mapper is, in percent, in plain words rather than "confidence"
    const sure = `the mapper is ${Math.round(Number(confidence) * 100)}% sure`;
    await expect.poll(async () => (await mappingState(page)).status, {message: 'how sure the mapper is',
      timeout: pollMs(180000)}).toContain(sure);
  }, {description: 'the Mapping section\'s status says how sure the mapper is: the confidence it gave (2 decimals) as a ' +
    'percent ("the mapper is 30% sure")'});

export const mappingWarn = Then('the Reaction panel\'s Mapping section shows a warn mentioning {string}',
  async (page: Page, word: string) => {
    await expect.poll(async () => (await mappingState(page)).flags.some((f) => f.level === 'warn' &&
      hasWords(f.text, [word])), {message: `a warn mentioning "${word}" in the Mapping section`,
    timeout: pollMs(180000)}).toBe(true);
  }, {description: 'the Mapping section shows a warn whose words hold the word (any case)'});

export const mappingWarn2 = Then('the Reaction panel\'s Mapping section shows a warn mentioning {string} and {string}',
  async (page: Page, a: string, b: string) => {
    await reactionPane(page);
    await expect.poll(async () => (await mappingState(page)).flags.some((f) => f.level === 'warn' &&
      hasWords(f.text, [a, b])), {message: `a warn mentioning "${a}" and "${b}" in the Mapping section`,
    timeout: pollMs(30000)}).toBe(true);
  }, {description: 'the Mapping section shows a warn whose words hold both words (any case)'});

export const mappingNoWarn = Then('the Reaction panel\'s Mapping section shows no warn', async (page: Page) => {
  await expect.poll(async () => mappedDone(await mappingState(page)), {message: 'the mapping finished',
    timeout: pollMs(180000)}).toBe(true);
  const s = await mappingState(page);
  expect(s.flags.filter((f) => f.level === 'warn'), 'no warn in the Mapping section').toEqual([]);
  // the status's words: how sure the mapper is, in percent, not "confidence"
  expect(s.status, 'the mapping\'s status (Crux mode: mapped, how sure the mapper is)').toMatch(/mapper is \d+% sure/);
}, {description: 'the mapping finished saying how sure the mapper is, and the Mapping section shows no warn'});

export const mappingNotShow = Then('the Reaction panel\'s Mapping section does not show {string}',
  async (page: Page, text: string) => {
    await expect.poll(async () => mappedDone(await mappingState(page)), {message: 'the mapping finished',
      timeout: pollMs(180000)}).toBe(true);
    expect((await mappingState(page)).status, 'the Mapping section').not.toContain(text);
  }, {description: 'once the mapping finished, its section says nothing with that text'});

export const mappingNote = Then('the Reaction panel\'s Mapping section shows a note mentioning {string}',
  async (page: Page, word: string) => {
    await expect.poll(async () => (await mappingState(page)).flags.some((f) => f.level === 'note' &&
      hasWords(f.text, [word])), {message: `a note mentioning "${word}" in the Mapping section`,
    timeout: pollMs(180000)}).toBe(true);
  }, {description: 'the Mapping section shows a note whose words hold the word (any case)'});

// ---------------------------------------------------------------------------------------------------- the engine-poison foot-gun

/* A getter of Crux Sketch throws after an edit when its engine has stopped (a trap, or it never loaded) or cannot write
   the drawing in the column's notation (the engine-poison foot-gun). No build stops the
   engine on demand, so the steps make the sketcher dialog's Crux getters throw in the page: own properties over the
   handle's four getters (smiles, molfile, molV3000, smarts), nothing else of it changed (no Proxy: its methods, destroy
   among them, keep their own `this`). Then a person's edit, as Crux's own 'user' change handler makes it: the host's
   explicit value (the cell's own SMILES) ends, and onChanged. */
async function crippleSketcherDialog(page: Page, error: {format: string, message: string}): Promise<void> {
  const find = `const ds = [...document.querySelectorAll('.d4-dialog')].filter((d) => d.querySelector('.crux-sketcher'));
    const root = ds[ds.length - 1]?.querySelector('[data-u2-name="cruxSketch"]');
    const impl = root ? DG.Widget.find(root) : null;`;
  await expect.poll(() => page.evaluate(new Function(`${find} return !!impl?.sketch && !impl.isRenderPending;`) as
    () => boolean), {message: 'the sketcher dialog\'s Crux ready', timeout: pollMs(30000)}).toBe(true);
  await page.evaluate(([body, err]) => new Function('err', body)(err), [`${find}
    const handle = impl.sketch;
    for (const n of ['smiles', 'molfile', 'molV3000', 'smarts'])
      Object.defineProperty(handle, n, {configurable: true, get: () => {
        throw Object.freeze({...err});
      }});
    // the host relays Crux's changes once its own init has resolved, a moment after Crux is ready (P6, 03:46: a change
    // told at once was lost): the edit is told again every 200 ms until the footer's slot answers (a block, or the
    // engine stop's warn), for at most 5 s
    const slot = ds[ds.length - 1].querySelector('.chem-footgun-slot');
    return (async () => {
      for (let i = 0; i < 25 && slot?.dataset.level !== 'block' && slot?.dataset.level !== 'warn'; i++) {
        impl.explicitMol = null;
        impl.onChanged.next(null);
        await new Promise((r) => setTimeout(r, 200));
      }
    })();`, error] as [string, {format: string, message: string}]);
}

export const engineStops = When('the Crux engine of the sketcher dialog stops', (page: Page) =>
  crippleSketcherDialog(page, {format: 'engine', message: 'Crux Sketch: the engine failed (RuntimeError: ' +
    'unreachable); reload the page to use the sketcher again'}),
{tier: 'api', description: 'in the page, the sketcher dialog\'s Crux getters throw crux-sketch\'s engine error ({format: \'engine\'}), as after a trap; then a change, as an edit tells it'});

export const enterInSketcherDialog = When('user presses Enter while the sketcher dialog has the focus',
  async (page: Page) => {
    const dialog = page.locator(dialogSelector('sketcher')).last();
    await dialog.locator('.d4-dialog-header, .d4-dialog-title').first().click();
    await page.keyboard.press('Enter');
    await page.waitForTimeout(700);
    const open = await page.locator(dialogSelector('sketcher')).count() > 0;
    console.log(`GUARD engine-poison: after Enter the sketcher dialog is ${open ? 'still open: the dialog did not ' +
      'take Enter while OK was disabled' : 'closed: the dialog took Enter, and the OK guard is what kept the cell'}`);
  }, {tier: 'ui', description: 'a click on the sketcher dialog\'s title (the focus in the dialog), then Enter; the log says whether the dialog took it'});

export const okDisabledOrNot = When('user clicks the OK button of the sketcher dialog, disabled or not',
  async (page: Page) => {
    const result = await page.evaluate((sel) => {
      const ds = [...document.querySelectorAll(sel)];
      const dialog = ds[ds.length - 1];
      if (!dialog)
        return 'not clicked: the dialog had closed';
      const ok = [...dialog.querySelectorAll('.d4-dialog-footer button')]
        .find((b) => /^\s*OK\s*$/i.test(b.textContent ?? '')) as HTMLButtonElement | undefined;
      if (!ok)
        return 'not clicked: no OK';
      const disabled = ok.disabled || ok.classList.contains('d4-disabled');
      ok.disabled = false;
      ok.classList.remove('d4-disabled', 'disabled');
      ok.click();
      return disabled ? 'clicked though disabled (enabled by hand)' : 'clicked';
    }, dialogSelector('sketcher'));
    console.log(`GUARD engine-poison: OK ${result}`);
    await page.waitForTimeout(700);
  }, {tier: 'ui', description: 'the sketcher dialog\'s OK clicked, enabled by hand first when disabled (what Enter or a stray handler could do); nothing when the dialog has closed'});

export const valueAsOpened = Then('the value in row {int} of {string} column is exactly what it was',
  async (page: Page, row: number, column: string) => {
    const was = opened.get(page)?.[column]?.[row - 1];
    expect(was, `the value row ${row} of "${column}" held as the table opened`).not.toBeUndefined();
    const now = await page.evaluate(([c, r]) => grok.shell.tv.dataFrame.col(c).get(r), [column, row - 1] as
      [string, number]);
    expect(now, `row ${row} of "${column}"`).toBe(was);
  }, {tier: 'api', description: 'the cell holds, character for character, the value the step that opened the table put there'});

