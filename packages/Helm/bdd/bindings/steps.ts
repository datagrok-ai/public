/* The steps only Helm defines: its readiness, a click on an empty spot of the editor canvas, typing
   into the editor's notation pane, and a sequence too long for the properties calculation. Everything
   generic — grid areas, the context menu, the context panel, dialogs, typing — is the library's
   (`grok-bdd list-steps`). */
import type {Page} from '@playwright/test';
import {Given, When} from '@datagrok-libraries/bdd';
import {ElementRef, expect, locate} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** Helm initializes on its first call (RDKit, the monomer library, the helper); the platform holds
 * the package's calls until then, so a service getter awaited here is that init made visible. */
export const helmInitialized = Given('the Helm package is initialized', async (page: Page) => {
  await page.evaluate(async () => { await grok.functions.call('Helm:getHelmHelper', {}); });
}, {tier: 'api', description: 'awaits Helm\'s init through the platform (its helper getter returns once initHelm has run)'});

/** Placing an armed monomer needs a click on the canvas away from every drawn monomer: the point
 * is picked from the monomers' own boxes, so a click on a monomer (which selects it) is never it. */
export const clickEmptyCanvas = When('user clicks on an empty spot of editor canvas', async (page: Page) => {
  const svg = await locate(page, 'editor canvas');
  const box = await svg.boundingBox();
  if (!box)
    throw new Error('editor canvas: not on the page');
  const atoms = await svg.locator('[data-testid^="canvas-atom-"]').evaluateAll((els) =>
    els.map((e) => e.getBoundingClientRect()).map((r) => ({x: r.x, y: r.y, w: r.width, h: r.height})));
  const clear = (x: number, y: number) => atoms.every((a) => x < a.x - 40 || x > a.x + a.w + 40 || y < a.y - 40 || y > a.y + a.h + 40);
  for (const [fx, fy] of [[0.85, 0.8], [0.15, 0.8], [0.85, 0.2], [0.15, 0.2], [0.5, 0.9]]) {
    const x = box.x + box.width * fx;
    const y = box.y + box.height * fy;
    if (clear(x, y)) {
      await page.mouse.click(x, y);
      return;
    }
  }
  throw new Error(`editor canvas: no spot 40 px clear of the ${atoms.length} drawn monomers`);
}, {tier: 'ui'});

/** The editor app takes Control+A for its own select-all (every monomer), so the library's typing,
 * which selects the old text with that key, appends to the notation instead of replacing it: the
 * text is selected the way a drag over it would, then typed key by key. */
export const replaceNotation = When('user replaces the text of {element} with {string}', async (page: Page, target: ElementRef, text: string) => {
  const loc = await locate(page, target);
  await loc.selectText();
  await page.keyboard.type(text);
  await expect(loc, `the text typed into ${target.phrase}`).toHaveText(text);
}, {tier: 'ui', description: 'for a contenteditable whose app owns Control+A; the text is typed, not committed'});

/** A peptide of N alanines: the properties calculation refuses a sequence over 1000 characters. */
export const setLongPeptide = When('user sets {string} column in row {int} to a peptide of {int} alanines', async (page: Page, column: string, row: number, n: number) => {
  await page.evaluate(([c, r, count]) => {
    grok.shell.t.col(c).set(r - 1, `PEPTIDE1{${Array(count).fill('A').join('.')}}$$$$`);
  }, [column, row, n] as [string, number, number]);
}, {tier: 'api'});
