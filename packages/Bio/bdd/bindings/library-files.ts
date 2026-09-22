/* The monomer library files as the package ships them, the standardized library the Match
   dialog's backend builds, and the readiness of the Manage Monomers view's sketcher. */
import type {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {callFunction, expect, pollMs, readResult} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const LIBRARIES = 'System:AppData/Bio/monomer-libraries/';

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
  const text: string = await page.evaluate((path) => grok.dapi.files.readAsText(path), LIBRARIES + name);
  await callFunction(page, 'Bio:standardiseMonomerLibrary', [['library', text]]);
}, {tier: 'api', description: 'Bio:standardiseMonomerLibrary over the shipped file; its JSON becomes the result the next steps read'});

export const resultHoldsMonomer = Then('the result should hold {string} monomer {string}', async (page: Page, polymerType: string, symbol: string) => {
  const found: string[] = await readResult(page, `((v) => Array.isArray(v) ? v : [])(typeof value === 'string' ? JSON.parse(value) : value).filter((m) => m.symbol === arg).map((m) => String(m.polymerType))`, symbol);
  expect(found, `polymer types of the "${symbol}" monomers in the result`).toContain(polymerType);
}, {description: 'a monomer list (a standardized library): some entry with that symbol and that polymer type'});

/** The Manage Monomers view hosts a monomer editor with a sketcher that mounts seconds after the
 * view (Ketcher: about ten on dev, over a minute on a stand serving a second worker). Closed
 * before that, the sketcher's late mount throws ResizeObserver TypeErrors into whatever runs
 * next, so a feature that opens the view lets the editor finish first — a Ketcher toolbar or,
 * for another backend, its canvas. */
export const monomerSketcherReady = Then('the monomer sketcher of the Manage Monomers view should be ready', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const root = document.querySelector('.monomer-manager-sketcher');
    if (!root)
      return 'no sketcher';
    const ketcher = root.querySelector('.Ketcher-root');
    return ketcher ? (ketcher.querySelectorAll('button').length > 5 ? 'ready' : 'mounting') : (root.querySelector('canvas') ? 'ready' : 'mounting');
  }), {timeout: pollMs(150000), message: 'the monomer editor sketcher'}).toBe('ready');
}, {description: 'the sketcher inside .monomer-manager-sketcher has its toolbar (Ketcher) or its canvas (other backends)'});
