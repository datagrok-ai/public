/* The generic step vocabulary: gestures (When) and outcomes (Then) over any element phrase.
   Expressions are cucumber expressions: `(on )` optional text, `in(to)` optional suffix, `be/become`
   alternation. Every definition is an exported const — the compiler imports it by name. */
import {readFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Given, Then, When} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {el} from '../../src/runtime/args.js';
import {expectCount, expectOptions, expectState, expectSwitched, expectText, expectValue, expectValueBetween, expectVisible, State} from '../../src/runtime/assertions.js';
import * as g from '../../src/runtime/gestures.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';
import {locate} from '../../src/runtime/locate.js';

export const clickOn = When('user clicks (on ){element}', (page: Page, target: ElementRef) => g.click(page, target), {tier: 'ui'});
export const doubleClickOn = When('user double-clicks (on ){element}', (page: Page, target: ElementRef) => g.dblclick(page, target), {tier: 'ui'});
export const rightClickOn = When('user right-clicks (on ){element}', (page: Page, target: ElementRef) => g.rightclick(page, target), {tier: 'ui'});
export const hoverOver = When('user hovers (over ){element}', (page: Page, target: ElementRef) => g.hover(page, target), {tier: 'ui'});
export const focusOn = When('user focuses (on ){element}', (page: Page, target: ElementRef) => g.focus(page, target), {tier: 'ui'});
export const typeInto = When('user types {string} in(to) {element}', (page: Page, text: string, target: ElementRef) => g.typeInto(page, target, text), {tier: 'ui'});
export const enterInto = When('user enters {string} in(to) {element}', (page: Page, text: string, target: ElementRef) => g.typeInto(page, target, text, true),
  {tier: 'ui', description: 'types and commits (Tab)'});
export const insertLine = When('user puts {string} on the first line of {element}',
  (page: Page, text: string, target: ElementRef) => g.insertLine(page, target, text),
  {tier: 'ui', description: 'a line typed into a code editor, whose document is not an input value'});
export const clearField = When('user clears {element}', (page: Page, target: ElementRef) => g.clear(page, target), {tier: 'ui'});
export const pressKey = When('user presses {key}', (page: Page, key: string) => g.press(page, key), {tier: 'ui'});
export const pressKeyIn = When('user presses {key} in {element}', (page: Page, key: string, target: ElementRef) => g.pressIn(page, target, key), {tier: 'ui'});
export const selectIn = When('user selects {string} in {element}', (page: Page, option: string, target: ElementRef) => g.select(page, target, option), {tier: 'ui'});
export const check = When('user checks {element}', (page: Page, target: ElementRef) => g.setChecked(page, target, true), {tier: 'ui'});
export const uncheck = When('user unchecks {element}', (page: Page, target: ElementRef) => g.setChecked(page, target, false), {tier: 'ui'});
export const toggle = When('user toggles {element}', (page: Page, target: ElementRef) => g.toggle(page, target), {tier: 'ui'});
export const switchOn = When('user switches on {element}', (page: Page, target: ElementRef) => g.setSwitched(page, target, true),
  {tier: 'ui', description: 'the switch that governs it — its own, or the one a parameter form puts beside it'});
export const switchOff = When('user switches off {element}', (page: Page, target: ElementRef) => g.setSwitched(page, target, false), {tier: 'ui'});
export const close = When('user closes {element}', (page: Page, target: ElementRef) => g.close(page, target), {tier: 'ui'});
export const expand = When('user expands {element}', (page: Page, target: ElementRef) => g.setExpanded(page, target, true),
  {tier: 'ui', description: 'tree nodes, accordion panes, dropdowns — anything with aria-expanded'});
export const collapse = When('user collapses {element}', (page: Page, target: ElementRef) => g.setExpanded(page, target, false), {tier: 'ui'});
export const dragSliderTo = When('user drags the slider of {element} to {float}',
  (page: Page, target: ElementRef, value: number) => g.dragSlider(page, target, value),
  {tier: 'ui', description: 'a real pointer drag along the track, to where the value lives on it'});
export const dragTo = When('user drags {element} to {element}', (page: Page, source: ElementRef, target: ElementRef) => g.drag(page, source, target), {tier: 'ui'});

/** `| element | value |` rows: selects choose an option, checkboxes take yes/no, everything else is typed. */
export const fillIn = When('user fills in:', async (page: Page, table: string[][]) => {
  for (const [phrase, value] of table) {
    const target = el(phrase);
    const editor = await g.editorOf(page, target);
    const shape = await editor.evaluate((e) => {
      const role = e.getAttribute('role') ?? e.querySelector('[role="switch"], [role="checkbox"], [role="radio"]')?.getAttribute('role');
      return `${e.tagName}:${(e as HTMLInputElement).type ?? ''}:${role ?? ''}`;
    }).catch(() => '');
    if (shape.startsWith('SELECT:'))
      await g.select(page, target, value);
    else if (/^INPUT:(checkbox|radio):|:(switch|checkbox|radio)$/.test(shape))
      await g.setChecked(page, target, /^(yes|true|on|checked|1)$/i.test(value));
    else
      await g.typeInto(page, target, value, true);
  }
}, {tier: 'ui'});

export const shouldBe = Then('{element} should be/become {state}', (page: Page, target: ElementRef, state: State) => expectState(page, target, state));
export const shouldNotBe = Then('{element} should not be/become {state}', (page: Page, target: ElementRef, state: State) => expectState(page, target, state, true));
export const shouldBecomeVisibleWithin = Then('{element} should become visible within {int} seconds',
  async (page: Page, target: ElementRef, seconds: number) => expectVisible(await locate(page, target), true, pollMs(seconds * 1000)),
  {description: 'for what a computation produces well past the usual budget (a search\'s hits): the budget is the scenario\'s claim about how long it may take'});
export const shouldContainText = Then('{element} should contain (the )text {string}', (page: Page, target: ElementRef, text: string) => expectText(page, target, text));
export const shouldNotContainText = Then('{element} should not contain (the )text {string}', (page: Page, target: ElementRef, text: string) => expectText(page, target, text, {negate: true}));
export const shouldHaveText = Then('{element} should have (the )text {string}', (page: Page, target: ElementRef, text: string) => expectText(page, target, text, {exact: true}));
export const shouldHaveValue = Then('{element} should have (the )value {string}', (page: Page, target: ElementRef, value: string) => expectValue(page, target, value));
export const shouldNotHaveValue = Then('{element} should not have (the )value {string}', (page: Page, target: ElementRef, value: string) => expectValue(page, target, value, true));
export const shouldBeSwitchedOn = Then('{element} should be switched on', (page: Page, target: ElementRef) => expectSwitched(page, target, true));
export const shouldBeSwitchedOff = Then('{element} should be switched off', (page: Page, target: ElementRef) => expectSwitched(page, target, false));
export const shouldHaveValueBetween = Then('{element} should have a value between {float} and {float}',
  (page: Page, target: ElementRef, lo: number, hi: number) => expectValueBetween(page, target, lo, hi),
  {description: 'a number a slider or a stepper arrives at, which no exact value would describe'});
export const visibleCount = Then('there should be {int} visible {element}', async (page: Page, count: number, target: ElementRef) =>
  expect((await locate(page, target)).filter({visible: true}), `visible ${target.phrase}`).toHaveCount(count),
{description: 'how many of the elements the phrase names are shown — one of a kind, never a duplicate'});

/** Visible counts remembered by phrase, for a claim against before (a search narrows a list whose
 * size depends on the stand: what matters is that it shrank). */
const rememberedVisible = new WeakMap<Page, Map<string, number>>();
const visibleNow = async (page: Page, target: ElementRef): Promise<number> => (await locate(page, target)).filter({visible: true}).count();
async function expectVisibleAgainstRemembered(page: Page, target: ElementRef, fewer: boolean): Promise<void> {
  const remembered = rememberedVisible.get(page)?.get(target.phrase);
  if (remembered === undefined)
    throw new Error(`no count of visible ${target.phrase} remembered: "user remembers the number of visible ${target.phrase}" first`);
  await expect.poll(async () => {
    const count = await visibleNow(page, target);
    return `${(fewer ? count < remembered : count > remembered) ? '' : 'not '}${fewer ? 'fewer' : 'more'} (${count} vs ${remembered})`;
  }, {message: `visible ${target.phrase} against the remembered count`}).toMatch(/^(fewer|more) \(/);
}
export const rememberVisibleCount = When('user remembers the number of visible {element}', async (page: Page, target: ElementRef) => {
  const counts = rememberedVisible.get(page) ?? new Map<string, number>();
  rememberedVisible.set(page, counts);
  counts.set(target.phrase, await visibleNow(page, target));
}, {description: 'how many are shown now, read once, for "fewer/more visible … than remembered" with the same phrase'});
export const visibleFewerThanRemembered = Then('there should be fewer visible {element} than remembered', (page: Page, target: ElementRef) =>
  expectVisibleAgainstRemembered(page, target, true), {description: 'strictly fewer than the remembered count — a list a search narrowed'});
export const visibleMoreThanRemembered = Then('there should be more visible {element} than remembered', (page: Page, target: ElementRef) =>
  expectVisibleAgainstRemembered(page, target, false), {description: 'strictly more than the remembered count — a list something added to'});
export const fillsParent = Then('{element} should fill its parent', async (page: Page, target: ElementRef) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  await expect.poll(() => loc.evaluate((e) => {
    const own = e.getBoundingClientRect();
    const host = e.parentElement!;
    const style = getComputedStyle(host);
    const width = host.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight);
    const height = host.clientHeight - parseFloat(style.paddingTop) - parseFloat(style.paddingBottom);
    return Math.abs(own.width - width) <= 1 && Math.abs(own.height - height) <= 1 ? 'fills' : `${Math.round(own.width)}×${Math.round(own.height)} in ${Math.round(width)}×${Math.round(height)}`;
  }), {message: `${target.phrase} against its parent`}).toBe('fills');
}, {description: 'as wide and as tall as the content box of the element it sits in, to a pixel'});
export const shouldOffer = Then('{element} should offer {string}', (page: Page, target: ElementRef, list: string) => expectOptions(page, target, list),
  {description: 'the choices of a dropdown, comma-separated, exactly and in this order; the blank option of a nullable dropdown is not a choice'});
export const shouldHaveItems = Then('{element} should have {int} item(s)', (page: Page, target: ElementRef, count: number) => expectCount(page, target, count));
export const shouldHaveRows = Then('{element} should have {int} row(s)', (page: Page, target: ElementRef, count: number) => expectCount(page, target, count));

/** `| element |` rows, every one checked for the state. */
export const followingShouldBe = Then('the following elements should be {state}:', async (page: Page, state: State, table: string[][]) => {
  for (const [phrase] of table)
    await expectState(page, el(phrase), state);
});

export const uploadThrough = When('user uploads {string} through {element}', (page: Page, file: string, target: ElementRef) => g.chooseFile(page, target, file),
  {tier: 'ui', description: 'clicks the element and answers the file chooser it opens with a file of the bdd project (a path under its root, "fixtures/lib.json")'});

export const clipboardContains = Then('the clipboard should contain (the )text {string}', async (page: Page, text: string) => {
  await expect.poll(() => g.readClipboard(page), {message: 'the clipboard text'}).toContain(text);
}, {description: 'what the page last copied (navigator.clipboard) — headless Chromium keeps a clipboard of its own'});

export const clipboardHas = Then('the clipboard should have (the )text {string}', async (page: Page, text: string) => {
  await expect.poll(() => g.readClipboard(page), {message: 'the clipboard text'}).toBe(text);
}, {description: 'exactly, whitespace included'});

const clipboardLines = async (page: Page): Promise<string[]> =>
  (await g.readClipboard(page)).split(/\r?\n/).filter((line) => line.length > 0);

export const clipboardLineCount = Then('the clipboard should hold {int} line(s)', async (page: Page, count: number) => {
  await expect.poll(async () => (await clipboardLines(page)).length, {message: 'the non-empty lines of the clipboard text'}).toBe(count);
}, {description: 'the non-empty lines of what the page last copied'});

export const clipboardLineValues = Then('line {int} of the clipboard should hold the tab-separated values {string}', async (page: Page, line: number, values: string) => {
  const want = values.split(/\s*,\s*/);
  await expect.poll(async () => ((await clipboardLines(page))[line - 1] ?? '').split('\t'),
    {message: `line ${line} of the clipboard, split at its tabs`}).toEqual(want);
}, {description: 'the values of one line of the clipboard (counted from 1, empty lines skipped), split at its tabs and compared with the comma-separated list, in order'});
export const pasteInto = When('user pastes {string} into {element}', async (page: Page, text: string, target: ElementRef) =>
  g.paste(page, await g.editorOf(page, target), text), {tier: 'ui', description: 'through the clipboard and the paste key over what the editor held; "\\n" is a line break'});

const rememberedClips = new WeakMap<Page, string[]>();

export const rememberClipboard = When('user remembers the clipboard text', async (page: Page) => {
  const text = await g.readClipboard(page);
  expect(text, 'the clipboard text to remember').not.toBe('');
  if (!rememberedClips.has(page))
    atFeatureEnd(page, async () => void rememberedClips.delete(page));
  rememberedClips.set(page, [...(rememberedClips.get(page) ?? []), text]);
}, {tier: 'api', description: 'kept with every text remembered before it in the feature'});

export const clipboardDiffers = Then('the clipboard text should differ from every remembered one', async (page: Page) => {
  const text = await g.readClipboard(page);
  const kept = rememberedClips.get(page) ?? [];
  expect(kept.length, 'texts remembered before').toBeGreaterThan(0);
  expect(kept.filter((k) => k === text).length, 'remembered texts equal to the clipboard').toBe(0);
});

export const clipboardImage = Then('the clipboard should hold a PNG image of at least {int} bytes', async (page: Page, bytes: number) => {
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  await expect.poll(() => page.evaluate(async () => {
    for (const item of await navigator.clipboard.read())
      if (item.types.includes('image/png'))
        return (await item.getType('image/png')).size;
    return 0;
  }), {message: 'the size of the PNG image on the clipboard'}).toBeGreaterThanOrEqual(bytes);
});

/** The state a scenario needs, rather than a gesture: `setExpanded` reads where the element is
 * first, so a group that is already open stays open — "user expands" on it would close it. */
export const isExpanded = Given('{element} is expanded', (page: Page, target: ElementRef) => g.setExpanded(page, target, true), {tier: 'ui'});

export const finishedUpdating = Then('{element} should have finished updating', async (page: Page, target: ElementRef) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  // ui.setUpdateIndicator shades the element it works in ("Updating..." with a loader) until the
  // work is done; a computation the click started (an MCS over a column) can take a while
  await expect.poll(() => loc.locator('.d4-update-shadow').count(),
    {timeout: pollMs(120000), message: `the update indicator ${target.phrase} shows`}).toBe(0);
}, {description: 'no "Updating..." shade over the element: what a click started in it has finished (up to two minutes)'});

// --- a file the page hands over, and handing it back ---------------------------------------------

/* The page saves a file (a form's "Save to file", an export, a cell action) through the browser's
   download. Every file the page downloads is kept in the temp directory, newest last: for the checks
   on its text and for the file chooser a later step answers with it. */
type Downloaded = {name: string; file: Promise<string>};
const downloads = new WeakMap<Page, Downloaded[]>();

/** The page's downloads; `fresh` starts the list afresh. The listener is attached once per page. */
function watched(page: Page, fresh: boolean): Downloaded[] {
  if (!downloads.has(page))
    page.on('download', (d) => {
      const file = join(tmpdir(), `bdd-${process.pid}-${Date.now()}-${d.suggestedFilename()}`);
      const saved = d.saveAs(file).then(() => file);
      // a download the page cancels rejects here, and only a step that reads that file reports it
      saved.catch(() => undefined);
      downloads.get(page)!.push({name: d.suggestedFilename(), file: saved});
    });
  if (fresh || !downloads.has(page))
    downloads.set(page, []);
  return downloads.get(page)!;
}

export const watchDownloads = Given('user watches downloads', async (page: Page) => {
  watched(page, true);
}, {tier: 'api', description: 'records the files the page downloads from then on, forgetting the ones before, so a file downloaded earlier cannot answer for one the scenario expects'});

export const downloadThrough = When('user downloads a file through {element}', async (page: Page, target: ElementRef) => {
  const before = watched(page, false).length;
  await g.click(page, target);
  await expect.poll(() => downloads.get(page)!.length, {timeout: 10000, message: `a download started by a click on ${target.phrase}`})
    .toBeGreaterThan(before);
}, {tier: 'ui', description: 'clicks the element and keeps the file the browser downloads, for "uploads the downloaded file"'});

async function lastDownloaded(page: Page): Promise<string> {
  const last = downloads.get(page)?.at(-1);
  if (!last)
    throw new Error('no file has been downloaded on this page');
  return last.file;
}

/** The browser numbers a file it has downloaded before ("smiles (2).sdf"), so the name a feature
 * gives matches those too. */
function named(name: string): (file: {name: string}) => boolean {
  const dot = name.lastIndexOf('.');
  const [stem, ext] = dot < 0 ? [name, ''] : [name.slice(0, dot), name.slice(dot)];
  const re = new RegExp(`^${stem.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}( \\(\\d+\\))?${ext.replace(/\./g, '\\.')}$`);
  return (file) => re.test(file.name);
}

async function textOf(page: Page, name: string): Promise<string> {
  let names = '';
  await expect.poll(() => {
    const list = downloads.get(page);
    if (!list)
      throw new Error('"user watches downloads" did not run in this scenario');
    names = list.map((d) => d.name).join(', ');
    return list.some(named(name));
  }, {message: `a download named "${name}"`}).toBe(true).catch(() => {
    throw new Error(`no download named "${name}"; downloaded: ${names || 'nothing'}`);
  });
  return readFileSync(await downloads.get(page)!.filter(named(name)).pop()!.file, 'utf8');
}

export const fileDownloaded = Then('a file {string} should have been downloaded', async (page: Page, name: string) => {
  await textOf(page, name);
});

export const downloadContains = Then('the downloaded file {string} should contain (the )text {string}', async (page: Page, name: string, text: string) => {
  expect(await textOf(page, name), `the text of "${name}"`).toContain(text);
});

const occurrences = async (page: Page, name: string, text: string): Promise<number> => (await textOf(page, name)).split(text).length - 1;

export const downloadCount = Then('the downloaded file {string} should contain {int} occurrences of {string}', async (page: Page, name: string, count: number, text: string) => {
  expect(await occurrences(page, name, text), `occurrences of "${text}" in "${name}"`).toBe(count);
}, {description: 'how many times the text appears in the file — a record terminator counts the records'});

export const downloadFewer = Then('the downloaded file {string} should contain fewer than {int} occurrences of {string}', async (page: Page, name: string, count: number, text: string) => {
  expect(await occurrences(page, name, text), `occurrences of "${text}" in "${name}"`).toBeLessThan(count);
});

export const downloadedContains = Then('the downloaded file should contain {string}', async (page: Page, text: string) => {
  const file = await lastDownloaded(page);
  expect(readFileSync(file, 'utf8'), `the downloaded file ${file}`).toContain(text);
}, {description: 'the file downloaded last'});

export const downloadedNotContains = Then('the downloaded file should not contain {string}', async (page: Page, text: string) => {
  const file = await lastDownloaded(page);
  expect(readFileSync(file, 'utf8'), `the downloaded file ${file}`).not.toContain(text);
}, {description: 'the file downloaded last'});

export const uploadDownloaded = When('user uploads the downloaded file through {element}', async (page: Page, target: ElementRef) =>
  g.chooseFile(page, target, await lastDownloaded(page)), {tier: 'ui', description: 'answers the file chooser the element opens with the file downloaded last'});

// --- code editors ---------------------------------------------------------------------------------

/* A CodeMirror document is not an input value: the text goes in at the caret, and it is read from
   the editor's own document — CM5 through the instance on its root, CM6 through the view the root
   carries — never from the rendered lines, which a long document virtualizes. */
async function codeOf(page: Page, target: ElementRef): Promise<string> {
  // a view keeps the editors of its other tabs in the page, hidden
  const loc = (await locate(page, target)).filter({visible: true}).first();
  return loc.evaluate((el: any) => {
    const root = el.classList?.contains('CodeMirror') || el.classList?.contains('cm-editor') ? el :
      el.querySelector('.CodeMirror, .cm-editor');
    if (root?.CodeMirror)
      return String(root.CodeMirror.getValue());
    const view = root?.cmView?.view ?? root?.querySelector?.('.cm-content')?.cmView?.view;
    if (view)
      return String(view.state.doc.toString());
    return String(root?.textContent ?? el.textContent ?? '');
  });
}

export const replaceCode = When('user replaces the code of {element} with {string}', async (page: Page, target: ElementRef, text: string) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  for (let attempt = 0; attempt < 2; attempt++) {
    await loc.click();
    await page.keyboard.press('ControlOrMeta+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type(text);
    if ((await codeOf(page, target)).trim() === text.trim())
      return;
  }
  expect((await codeOf(page, target)).trim(), `the code of ${target.phrase} after it was replaced`).toBe(text.trim());
}, {tier: 'ui', description: 'select-all and type at the caret, the document read back (retyped once when an editor mounting late ate keys)'});

export const holdsCode = Then('{element} should hold the code {string}', async (page: Page, target: ElementRef, text: string) => {
  await expect.poll(async () => (await codeOf(page, target)).trim(), {message: `the code of ${target.phrase}`}).toBe(text.trim());
}, {description: 'the editor document, exactly (trimmed)'});

export const appendToEditor = When('user appends {string} to {element}', async (page: Page, text: string, target: ElementRef) => {
  const editor = (await locate(page, target)).first();
  await editor.click();
  await page.keyboard.press('ControlOrMeta+End');
  await page.keyboard.press('Enter');
  await page.keyboard.type(text);
  await expect(editor, `${target.phrase} after the text was appended`).toContainText(text);
}, {tier: 'ui', description: 'a new last line typed into a code editor, so its dirty flag sees a real edit'});

/** A Dart text area (`ui.textInput` multiline) is a bare `<textarea class="ui-input-editor">` that
 * no input kind reaches, and its text is its value. */
export const textAreaHolds = Then('the text area of {element} should hold {string}', async (page: Page, target: ElementRef, text: string) => {
  const area = (await locate(page, target)).first().locator('textarea').first();
  await expect(area, `the text area of ${target.phrase}`).toHaveValue(text, {timeout: 15000});
}, {description: 'the value of the first <textarea> inside the element'});

// --- secrets, local files and the browser's own dialogs --------------------------------------------

const exposed = new WeakSet<Page>();

/** An environment variable (a password, an API key) put into the element's editor without passing
 * through the test's own calls: the page asks a function the test exposed, so neither the step, its
 * failure, the trace nor a guide carries the value. */
export const enterSecret = When('user enters the {word} secret into {element}', async (page: Page, variable: string, target: ElementRef) => {
  if (!process.env[variable])
    throw new Error(`the ${variable} environment variable is not set — a @needs-credentials scenario needs it`);
  if (!exposed.has(page)) {
    await page.exposeFunction('__bddSecret', (name: string) => process.env[name] ?? null);
    exposed.add(page);
  }
  const editor = await g.editorOf(page, target);
  const entered = await editor.evaluate(async (input, name) => {
    const value = await (window as any).__bddSecret(name);
    const field = input as HTMLInputElement;
    field.focus();
    field.value = value ?? '';
    field.dispatchEvent(new Event('input', {bubbles: true}));
    field.dispatchEvent(new Event('change', {bubbles: true}));
    field.blur();
    return field.value.length > 0;
  }, variable);
  expect(entered, `${variable} entered into ${target.phrase}`).toBe(true);
}, {tier: 'ui', description: 'the value is never printed — in a failure, a trace or a guide; the step fails naming the variable when it is unset'});

export const openLocalFile = When('user opens the local file {string}', async (page: Page, file: string) => {
  const chooser = page.waitForEvent('filechooser', {timeout: 10000});
  await page.keyboard.press('ControlOrMeta+O');
  await (await chooser).setFiles(join(process.env.BDD_ROOT ?? process.cwd(), file));
}, {tier: 'ui', description: 'Ctrl/Cmd+O and the file chooser it opens, answered with a file of the bdd project — the keyboard way in to what a drop from the desktop does'});

const alerts = new WeakMap<Page, string[]>();

export const recordAlerts = Given('browser alerts are recorded', async (page: Page) => {
  if (alerts.has(page))
    return;
  const list: string[] = [];
  alerts.set(page, list);
  page.on('dialog', async (dialog) => {
    list.push(dialog.message());
    await dialog.dismiss().catch(() => undefined);
  });
}, {tier: 'api', description: 'every native alert the page raises from now on is recorded and dismissed (an undismissed one blocks every later step)'});

export const alertShown = Then('the browser should have shown the alert {string}', async (page: Page, text: string) => {
  await expect.poll(() => alerts.get(page) ?? [], {message: 'the native alerts the page raised'}).toContain(text);
});

export const computedStyleContains = Then('the {string} style of {element} should contain {string}', async (page: Page, property: string, target: ElementRef, text: string) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  await expect.poll(() => loc.evaluate((e, p) => getComputedStyle(e).getPropertyValue(p), property),
    {message: `the "${property}" style of ${target.phrase}`}).toContain(text);
}, {description: 'the computed CSS value the browser rendered the element with (font-family, font-size, …), by substring'});
