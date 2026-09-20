/* The generic step vocabulary: gestures (When) and outcomes (Then) over any element phrase.
   Expressions are cucumber expressions: `(on )` optional text, `in(to)` optional suffix, `be/become`
   alternation. Every definition is an exported const — the compiler imports it by name. */
import {readFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
import {Given, Then, When} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {el} from '../../src/runtime/args.js';
import {expectCount, expectOptions, expectState, expectSwitched, expectText, expectValue, expectValueBetween, State} from '../../src/runtime/assertions.js';
import * as g from '../../src/runtime/gestures.js';
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
  {description: 'the choices of a dropdown, comma-separated, exactly and in this order'});
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

export const pasteInto = When('user pastes {string} into {element}', async (page: Page, text: string, target: ElementRef) =>
  g.paste(page, await g.editorOf(page, target), text), {tier: 'ui', description: 'through the clipboard and the paste key over what the editor held; "\\n" is a line break'});

/** The state a scenario needs, rather than a gesture: `setExpanded` reads where the element is
 * first, so a group that is already open stays open — "user expands" on it would close it. */
export const isExpanded = Given('{element} is expanded', (page: Page, target: ElementRef) => g.setExpanded(page, target, true), {tier: 'ui'});

// --- a file the page hands over, and handing it back ---------------------------------------------

/* The page saves a file (a form's "Save to file", an export) through the browser's download; the
   step keeps it in the temp directory for the file chooser a later step answers with it. One file
   per worker: the next download replaces it. */
let downloaded: string | undefined;

function downloadedFile(): string {
  if (!downloaded)
    throw new Error('no file has been downloaded in this worker');
  return downloaded;
}

export const downloadThrough = When('user downloads a file through {element}', async (page: Page, target: ElementRef) => {
  const download = page.waitForEvent('download', {timeout: 10000});
  await g.click(page, target);
  const file = await download;
  downloaded = join(tmpdir(), `bdd-${process.pid}-${Date.now()}-${file.suggestedFilename()}`);
  await file.saveAs(downloaded);
}, {tier: 'ui', description: 'clicks the element and keeps the file the browser downloads, for "uploads the downloaded file"'});

export const downloadedContains = Then('the downloaded file should contain {string}', async (page: Page, text: string) => {
  expect(readFileSync(downloadedFile(), 'utf8'), `the downloaded file ${downloaded}`).toContain(text);
});

export const downloadedNotContains = Then('the downloaded file should not contain {string}', async (page: Page, text: string) => {
  expect(readFileSync(downloadedFile(), 'utf8'), `the downloaded file ${downloaded}`).not.toContain(text);
});

export const uploadDownloaded = When('user uploads the downloaded file through {element}', (page: Page, target: ElementRef) =>
  g.chooseFile(page, target, downloadedFile()), {tier: 'ui', description: 'answers the file chooser the element opens with the file the last "downloads a file" kept'});

export const computedStyleContains = Then('the {string} style of {element} should contain {string}', async (page: Page, property: string, target: ElementRef, text: string) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  await expect.poll(() => loc.evaluate((e, p) => getComputedStyle(e).getPropertyValue(p), property),
    {message: `the "${property}" style of ${target.phrase}`}).toContain(text);
}, {description: 'the computed CSS value the browser rendered the element with (font-family, font-size, …), by substring'});
