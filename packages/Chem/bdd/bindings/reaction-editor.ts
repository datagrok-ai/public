// The reaction cell editor (src/utils/reactions/reaction-cell-editor.ts): its dialog, the arrows Crux draws in it, and
// whether a reaction cell still holds a reaction after OK. Crux's box in the dialog is a widget named cruxSketch, so its
// atoms and bonds are areas of the library's "crux sketcher widget"; "should differ from what was remembered" is beside
// its byte-for-byte twin in crux.ts. What the editor writes for the corpus and Chem's templates is Chem's package tests
// ('Reactions: editor'). The demo clips' steps follow: the gear's toolset items, the Reactions app's text box and
// drawing, and a reaction's reactants and products read by RDKit (a cell's, the app's text).
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

element('reaction editor dialog', {selector: '.d4-dialog:has(.chem-reaction-cell-editor)',
  description: 'the reaction cell editor: the dialog a double-click on a reaction cell opens (Crux mode: Crux\'s sketcher on ' +
    'the cell\'s reaction; RDKit mode: its text)'});

export const editorArrows = Then('the reaction editor should draw {int} reaction arrow(s)', async (page: Page, n: number) => {
  await expect.poll(() => page.evaluate(() => {
    const root = document.querySelector('.chem-reaction-cell-editor');
    const status = root ? DG.Widget.find(root)?.getWidgetStatus() : null;
    return status?.values?.ready ? status.values.arrows : null;
  }), {message: 'the arrows the reaction editor draws', timeout: pollMs(15000)}).toBe(n);
}, {tier: 'api', description: 'the reaction arrows Crux draws in the reaction cell editor (its status\'s `arrows` reading)'});

/** A row out of sight brought into view (three wheel notches over cell 4 once put row 10 past the top, and the
 * double-click on its cell timed out): scrolled to as a person drags the scroll bar to it. */
export const scrollGridToRow = When('user scrolls the grid to row {int} of the {string} column',
  async (page: Page, row: number, column: string) => {
    await page.evaluate(([r, c]) => {
      const grid = grok.shell.tv?.grid;
      if (!grid?.dataFrame.col(c))
        throw new Error(`the grid of the current view has no "${c}" column`);
      grid.scrollToCell(c, r - 1);
    }, [row, column] as [number, string]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'grid.scrollToCell on the row (from 1) of the column, then the viewers settle'});

export const editorArranged = Then(
  'the reaction editor should draw the reaction in at least {int} rows, its bonds at least {int} px long',
  async (page: Page, rows: number, px: number) => {
    await expect.poll(() => page.evaluate(([r, p]) => {
      const root = document.querySelector('.chem-reaction-cell-editor');
      const status = root ? DG.Widget.find(root)?.getWidgetStatus() : null;
      return status?.values?.ready ? status.values.rows >= r && status.values.bondLength >= p : null;
    }, [rows, px] as [number, number]), {message: 'the reaction editor\'s drawing arranged: its rows and bond length',
      timeout: pollMs(15000)}).toBe(true);
  }, {tier: 'api', description: 'the rows the reaction editor\'s arrows stand in and its drawn bond length (its status\'s ' +
    '`rows` and `bondLength` readings): a route opened arranged (decision 12)'});

export const readsAsReaction = Then('Chem should read the value of {string} column in row {int} as a reaction',
  async (page: Page, column: string, row: number) => {
    const read = await page.evaluate(async ([c, r]) => {
      const text = String(grok.shell.t.col(c).get(r - 1) ?? '');
      const rdkit = await grok.functions.call('Chem:getRdKitModule');
      let rxn: any = null;
      try {
        rxn = rdkit.get_rxn(text);
        return rxn !== null; // a molecule's text has no arrow: get_rxn refuses it
      } catch {
        return false;
      } finally {
        rxn?.delete();
      }
    }, [column, row] as [string, number]);
    expect(read, `${column} in row ${row} read as a reaction`).toBe(true);
  }, {tier: 'api', description: 'the cell is reaction text (RDKit\'s get_rxn reads it): not a molecule, as RXN-16 would have it'});

// ---------------------------------------------------------------- the demo clips (features/demos/reactions-edited-crux.feature)

// The gear's toolset items are gone (crux-sketch d0b8bb9): the host sets the toolset, the reaction toolset's bar has N.
element('reactions app text', {selector: '.chem-reactions-app-text',
  description: 'the Reactions app\'s text box, under "Text" below its sketcher: the reaction\'s text, which the sketcher follows'});
element('reactions app drawing', {selector: '.chem-reactions-app-drawing',
  description: 'the Reactions app\'s right pane: its title and the reaction drawn for publication'});

export const openReactionsApp = Given('user opens Chem\'s Reactions app', async (page: Page) => {
  await page.evaluate(async () => {
    const view = await grok.functions.call('Chem:reactionsApp');
    grok.shell.addView(view);
  });
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Reaction Sketcher');
}, {tier: 'api', description: 'runs Chem\'s Reactions app (its function `reactionsApp`, meta.role app, which the library\'s ' +
  '"user opens the {string} app" does not find: it looks for the app tag) and shows its view; done when that view is current'});

/** In the page: the reactants and the products of one reaction step's text (a CX part left out), each component as
 * RDKit writes it canonically, sorted; agents left out. Null for text that is not one step. */
function sidesOf(page: Page, text: string): Promise<{reactants: string[], products: string[]} | null> {
  return page.evaluate(async (t) => {
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const canonical = (smiles: string) => {
      let mol: any = null;
      try {
        mol = rdkit.get_mol(smiles);
        return mol?.is_valid() ? mol.get_smiles() : `(unread) ${smiles}`;
      } catch {
        return `(unread) ${smiles}`;
      } finally {
        mol?.delete();
      }
    };
    const parts = t.split(' |')[0].trim().split('>');
    if (parts.length !== 3)
      return null;
    const side = (part: string) => part === '' ? [] : part.split('.').map(canonical).sort();
    return {reactants: side(parts[0]), products: side(parts[2])};
  }, text);
}

/** Polls `read` (a reaction's text) until its reactants and products are the ones named, by RDKit. */
async function expectSides(page: Page, read: () => Promise<string>, reactants: string, products: string, what: string) {
  const want = JSON.stringify(await sidesOf(page, `${reactants}>>${products}`));
  let text = '';
  try {
    await expect.poll(async () => JSON.stringify(await sidesOf(page, text = await read())),
      {message: `${what}: the reactants and the products`, timeout: pollMs(10000)}).toBe(want);
  } catch (e) {
    throw new Error(`${what}: the text read is ${text}\n${e instanceof Error ? e.message : e}`);
  }
}

export const cellSides = Then('the reaction in row {int} of {string} column should have the reactants {string} and the products {string}',
  (page: Page, row: number, column: string, reactants: string, products: string) =>
    expectSides(page, () => page.evaluate(([c, r]) => String(grok.shell.t.col(c).get(r - 1) ?? ''),
      [column, row] as [string, number]), reactants, products, `${column} in row ${row}`),
  {tier: 'api', description: 'the cell\'s reaction (one step) has these reactants and products, each read by RDKit in any SMILES, the agents left out'});

export const appSides = Then('the Reactions app\'s text should have the reactants {string} and the products {string}',
  (page: Page, reactants: string, products: string) =>
    expectSides(page, () => page.evaluate(() =>
      (document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null)?.value ?? ''),
    reactants, products, 'the Reactions app\'s text'),
  {tier: 'api', description: 'the reaction in the Reactions app\'s text box (one step) has these reactants and products, each read by RDKit in any SMILES, the agents left out'});

/* For crux.ts, after `cellAsRemembered` (the cells it remembers are its own):

export const cellDiffers = Then('the value of {string} column in row {int} should differ from what was remembered',
  async (page: Page, column: string, row: number) => {
    const before = cells.get(page);
    if (before === undefined)
      throw new Error('no cell remembered: "user remembers the value of … column in row …" first');
    await expect.poll(() => page.evaluate(([c, r]) => String(grok.shell.t.col(c).get(r - 1) ?? ''), [column, row] as [string, number]),
      {message: `${column} in row ${row}`}).not.toBe(before);
  }, {description: 'the cell\'s text is no longer the string remembered'});
*/
