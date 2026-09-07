/* The generic step vocabulary: gestures (When) and outcomes (Then) over any element phrase.
   Expressions are cucumber expressions: `(on )` optional text, `in(to)` optional suffix, `be/become`
   alternation. Every definition is an exported const — the compiler imports it by name. */
import {expect, Page} from '@playwright/test';
import {Given, Then, When} from '../../src/registry.js';
import type {ElementRef} from '../../src/runtime/args.js';
import {el} from '../../src/runtime/args.js';
import {expectCount, expectState, expectText, expectValue, State} from '../../src/runtime/assertions.js';
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
export const clearField = When('user clears {element}', (page: Page, target: ElementRef) => g.clear(page, target), {tier: 'ui'});
export const pressKey = When('user presses {key}', (page: Page, key: string) => g.press(page, key), {tier: 'ui'});
export const pressKeyIn = When('user presses {key} in {element}', (page: Page, key: string, target: ElementRef) => g.pressIn(page, target, key), {tier: 'ui'});
export const selectIn = When('user selects {string} in {element}', (page: Page, option: string, target: ElementRef) => g.select(page, target, option), {tier: 'ui'});
export const check = When('user checks {element}', (page: Page, target: ElementRef) => g.setChecked(page, target, true), {tier: 'ui'});
export const uncheck = When('user unchecks {element}', (page: Page, target: ElementRef) => g.setChecked(page, target, false), {tier: 'ui'});
export const toggle = When('user toggles {element}', (page: Page, target: ElementRef) => g.toggle(page, target), {tier: 'ui'});
export const open = When('user opens {element}', (page: Page, target: ElementRef) => g.click(page, target), {tier: 'ui'});
export const close = When('user closes {element}', (page: Page, target: ElementRef) => g.close(page, target), {tier: 'ui'});
export const expand = When('user expands {element}', (page: Page, target: ElementRef) => g.setExpanded(page, target, true),
  {tier: 'ui', description: 'tree nodes, accordion panes, dropdowns — anything with aria-expanded'});
export const collapse = When('user collapses {element}', (page: Page, target: ElementRef) => g.setExpanded(page, target, false), {tier: 'ui'});
export const dragTo = When('user drags {element} to {element}', (page: Page, source: ElementRef, target: ElementRef) => g.drag(page, source, target), {tier: 'ui'});
export const scrollTo = When('user scrolls to {element}', (page: Page, target: ElementRef) => g.scrollTo(page, target), {tier: 'ui'});
export const navigateTo = When('user navigates to {string}', (page: Page, url: string) => page.goto(url, {waitUntil: 'domcontentloaded'}).then(() => undefined),
  {tier: 'ui', description: 'a URL or a path on the stand'});
export const reloadPage = When('user reloads the page', (page: Page) => page.reload({waitUntil: 'domcontentloaded'}).then(() => undefined), {tier: 'ui'});

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
export const shouldHaveItems = Then('{element} should have {int} item(s)', (page: Page, target: ElementRef, count: number) => expectCount(page, target, count));
export const shouldHaveRows = Then('{element} should have {int} row(s)', (page: Page, target: ElementRef, count: number) => expectCount(page, target, count));
export const shouldHaveTabs = Then('{element} should have {int} tab(s)', (page: Page, target: ElementRef, count: number) => expectCount(page, target, count));

/** `| element |` rows, every one checked for the state. */
export const followingShouldBe = Then('the following elements should be {state}:', async (page: Page, state: State, table: string[][]) => {
  for (const [phrase] of table)
    await expectState(page, el(phrase), state);
});

export const waitFor = Given('user waits for {element}', async (page: Page, target: ElementRef) => {
  await expect(await locate(page, target)).toBeVisible();
}, {description: 'until visible — the only wait the vocabulary has'});
