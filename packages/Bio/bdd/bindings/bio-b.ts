/* Readings the gap round needed: a library file as the package ships it, the standardized library
   the Match dialog's backend builds, and a numbering run's position names. */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {readResult} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const LIBRARIES = 'System:AppData/Bio/monomer-libraries/';

/** The HELM library schema: a JSON array of monomers, each with a symbol and a structure. */
export const libraryFileConforms = Then('the {string} monomer library file should list monomers with a symbol and a structure each',
  async (page: Page, name: string) => {
    const facts: {count: number, bad: string[]} = await page.evaluate(async (path) => {
      const monomers = JSON.parse(await grok.dapi.files.readAsText(path));
      if (!Array.isArray(monomers))
        throw new Error(`${path} is not a JSON array of monomers`);
      const bad = monomers.filter((m: any) => !m?.symbol || !(m.molfile || m.smiles)).map((m: any) => String(m?.symbol ?? JSON.stringify(m).slice(0, 40)));
      return {count: monomers.length, bad};
    }, LIBRARIES + name);
    expect(facts.count, `monomers in ${name}`).toBeGreaterThan(0);
    expect(facts.bad, `monomers of ${name} without a symbol or a molfile/smiles`).toEqual([]);
  }, {tier: 'api', description: 'the file under System:AppData/Bio/monomer-libraries: a non-empty array, every monomer with a symbol and a molfile or smiles'});

export const standardiseLibrary = When('user standardises the {string} monomer library', async (page: Page, name: string) => {
  await page.evaluate(async (path) => {
    const text = await grok.dapi.files.readAsText(path);
    const result = await grok.functions.call('Bio:standardiseMonomerLibrary', {library: text});
    (window as any).__bddLastResult = {name: 'Bio:standardiseMonomerLibrary', value: JSON.parse(result)};
  }, LIBRARIES + name);
}, {tier: 'api', description: 'Bio:standardiseMonomerLibrary over the shipped file; the parsed JSON becomes the result the next steps read'});

export const resultHoldsMonomer = Then('the result should hold {string} monomer {string}', async (page: Page, polymerType: string, symbol: string) => {
  const found: string[] = await readResult(page, `(Array.isArray(value) ? value : []).filter((m) => m.symbol === arg).map((m) => String(m.polymerType))`, symbol);
  expect(found, `polymer types of the "${symbol}" monomers in the result`).toContain(polymerType);
}, {description: 'a monomer list (a standardized library): some entry with that symbol and that polymer type'});

export const positionNamesCount = Then('{string} column should list at least {int} position names', async (page: Page, column: string, count: number) => {
  const names: string[] = await page.evaluate((c) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${grok.shell.t.name}`);
    return String(col.getTag('.positionNames') ?? '').split(',').map((s: string) => s.trim()).filter((s: string) => s !== '');
  }, column);
  expect(names.length, `position names in the .positionNames tag of "${column}"`).toBeGreaterThanOrEqual(count);
}, {description: 'the comma-separated .positionNames tag a numbering run writes on the aligned column'});

/** The Manage Monomers view hosts a monomer editor with a sketcher that mounts seconds after the
 * view (Ketcher: about ten on dev). Closed before that, the sketcher's late mount throws
 * ResizeObserver TypeErrors into whatever runs next, so a feature that opens the view lets the
 * editor finish first — a Ketcher toolbar or, for another backend, its canvas. */
export const monomerSketcherReady = Then('the monomer sketcher of the Manage Monomers view should be ready', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const root = document.querySelector('.monomer-manager-sketcher');
    if (!root)
      return 'no sketcher';
    const ketcher = root.querySelector('.Ketcher-root');
    return ketcher ? (ketcher.querySelectorAll('button').length > 5 ? 'ready' : 'mounting') : (root.querySelector('canvas') ? 'ready' : 'mounting');
  }), {timeout: 60000, message: 'the monomer editor sketcher'}).toBe('ready');
}, {description: 'the sketcher inside .monomer-manager-sketcher has its toolbar (Ketcher) or its canvas (other backends)'});
