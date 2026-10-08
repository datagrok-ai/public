/* Crux Sketch, Chem's own molecule sketcher (src/crux/crux-sketcher.ts), as Chem's features drive it. The sketcher's
   widget ("crux sketcher widget": its "atom N" and "bond N" areas, its readings), its controls by test id and the
   molecule readings are the library's `molecules` tier (bdd.config.json), shared with the other packages whose hosts
   open Crux. Here is what is Chem's own: the drawing a sketcher shows in place of itself (a filter card's, a pane's),
   Ketcher's canvas after a switch, and what the platform does with the sketcher that no gesture shows: the change
   events of the sketcher a step opens next, counted from its creation, and the molecules made the current object. */
import {Page} from '@playwright/test';
import {element, Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, el, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

element('sketcher thumbnail', {selector: '.chem-external-sketcher-canvas',
  description: 'the drawing a sketcher shows in place of itself until clicked (a filter card, a pane\'s scaffold): its Clear button shows while the pointer is over it'});

element('Ketcher canvas', {selector: '.d4-dialog .Ketcher-root [data-testid="ketcher-canvas"]',
  description: 'the canvas of Ketcher in a sketcher dialog, once Ketcher has loaded into it (a switch away from Crux lands there)'});

// ---------------------------------------------------------------- the next sketcher, from its creation

declare global {
  interface Window {
    __cruxNext?: {changes: number, created: boolean, cleared: boolean | null};
    __cruxCurrent?: {semType: string | null, value: string}[];
  }
}

/** In the page: watches for the next Crux sketcher to enter it, and calls `found(impl, host)` in the task it does,
 * before its init resolves (the host puts the sketcher's root in its box, then awaits its init). */
function watchNextCrux(page: Page, clearWhenReady: boolean): Promise<void> {
  return page.evaluate((clear) => {
    window.__cruxNext = {changes: 0, created: false, cleared: clear ? false : null};
    const observer = new MutationObserver((records) => {
      for (const r of records) {
        for (const n of Array.from(r.addedNodes)) {
          if (!(n instanceof HTMLElement) || !n.classList.contains('crux-sketcher'))
            continue;
          const impl = DG.Widget.find(n);
          if (!impl)
            continue;
          observer.disconnect();
          window.__cruxNext!.created = true;
          impl.onChanged.subscribe(() => window.__cruxNext!.changes++);
          if (!clear)
            return;
          // the sketcher's host (DG.chem.Sketcher), the nearest widget round the sketcher's root
          let host: any = null;
          for (let e = n.parentElement; e !== null && host === null; e = e.parentElement) {
            const w = DG.Widget.find(e);
            if (w && w.onSketcherReady)
              host = w;
          }
          const sub = host.onSketcherReady.subscribe((ready: any) => {
            if (ready !== impl)
              return;
            sub.unsubscribe();
            // the task after the announcement: after anything the host does on it in its own microtasks
            setTimeout(() => {
              (n.querySelector('crux-sketch')!.shadowRoot!.querySelector('[data-testid="toolbar.clear"]') as HTMLElement).click();
              window.__cruxNext!.cleared = true;
            }, 0);
          });
          return;
        }
      }
    });
    observer.observe(document.body, {childList: true, subtree: true});
  }, clearWhenReady);
}

export const countNextCrux = Given('the change events of the next Crux sketcher are counted', (page: Page) =>
  watchNextCrux(page, false),
{tier: 'api', description: 'the next Crux sketcher to enter the page has its onChanged counted from its creation, before its init resolves: one shown with a molecule has fired once'});

export const clearNextCrux = Given('the next Crux sketcher clears its canvas the task after its host says it is ready', (page: Page) =>
  watchNextCrux(page, true),
{tier: 'api', description: 'a user\'s edit at the earliest moment: Clear canvas pressed in the task after the host\'s onSketcherReady for it, after whatever the host does on that announcement in its own microtasks'});

export const cruxChanges = Then('the Crux sketcher should have fired {int} change event(s)', async (page: Page, n: number) => {
  const read = () => page.evaluate(() => window.__cruxNext ?? null);
  const first = await read();
  expect(first?.created, 'no Crux sketcher entered the page since the count began').toBe(true);
  // a count is final once the sketcher has nothing on its way
  await viewers.settle(page, el('crux sketcher widget')).catch(() => undefined);
  await expect.poll(async () => (await read())?.changes, {message: 'the Crux sketcher\'s change events'}).toBe(n);
}, {description: 'onChanged of the sketcher the count follows, read once it has nothing on its way'});

export const cruxCleared = Then('the Crux sketcher should have cleared its canvas', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => window.__cruxNext?.cleared ?? null), {message: 'Clear canvas pressed'}).toBe(true);
});

// ---------------------------------------------------------------- the molblock the sketcher writes

/** The atoms' coordinate lines of the host's molblock (DG.chem.Sketcher.getMolFile of the Crux in the page). */
function hostCoordinates(page: Page): Promise<string> {
  return page.evaluate(() => {
    const root = Array.from(document.querySelectorAll('.crux-sketcher')).pop();
    const impl = root ? DG.Widget.find(root) : null;
    if (!impl)
      throw new Error('no Crux sketcher in the page');
    const molFile: string = impl.molFile;
    return molFile.split('\n').filter((l) => /^\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+ [A-Z*]/.test(l))
      .map((l) => l.slice(0, 30)).join('|');
  });
}

const remembered = new WeakMap<Page, string>();

export const rememberCoordinates = When('user remembers the coordinates of the Crux sketcher\'s molblock', async (page: Page) => {
  remembered.set(page, await hostCoordinates(page));
}, {tier: 'api', description: 'the x, y, z of every atom of the molFile the sketcher writes'});

export const newCoordinates = Then('the Crux sketcher\'s molblock should have new coordinates', async (page: Page) => {
  const before = remembered.get(page);
  if (before === undefined)
    throw new Error('no coordinates remembered: "user remembers the coordinates of the Crux sketcher\'s molblock" first');
  await expect.poll(() => hostCoordinates(page), {message: 'the atoms\' coordinates in the sketcher\'s molFile'}).not.toBe(before);
}, {description: 'an atom\'s x, y or z in the molFile differs from the remembered one'});

// ---------------------------------------------------------------- a gesture that ends on a dialog's button

export const dragThenClick = When('user drags the {string} area of {widget} by {int} pixels to the {word} and at once clicks on {element}',
  async (page: Page, area: string, target: ElementRef, px: number, direction: string, button: ElementRef) => {
    const at = viewers.centerOf(await viewers.hitArea(page, target, area, true));
    const {dx, dy} = viewers.dragDelta(px, direction);
    const box = await (await locate(page, button)).filter({visible: true}).first().boundingBox();
    if (box === null)
      throw new Error(`${button.phrase} takes no space`);
    await page.mouse.move(at.x, at.y);
    await page.mouse.down();
    await page.mouse.move(at.x + dx / 2, at.y + dy / 2, {steps: 4});
    await page.mouse.move(at.x + dx, at.y + dy, {steps: 4});
    await page.mouse.up();
    // no wait and no read between the release and the press
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  }, {tier: 'ui', description: 'the drag of an area, and a click on the element right after the release (its place read before the drag)'});

// ---------------------------------------------------------------- the current object

export const watchCurrent = Given('the molecules made the current object are counted', (page: Page) =>
  page.evaluate(() => {
    window.__cruxCurrent = [];
    grok.events.onCurrentObjectChanged.subscribe(() => {
      const o = grok.shell.o;
      if (window.__cruxCurrent && o?.semType === DG.SEMTYPE.MOLECULE && typeof o.value === 'string' && o.value.trim() !== '')
        window.__cruxCurrent.push({semType: o.semType, value: o.value});
    });
  }),
{tier: 'api', description: 'from now on, each change of the current object to a Molecule semantic value with a molecule in it (what a sketcher\'s host sets on an edit); a dialog the platform makes current, or an empty cell, is not one'});

export const currentCount = Then('a molecule should have become the current object {int} time(s)', async (page: Page, n: number) => {
  await viewers.settle(page, el('crux sketcher widget')).catch(() => undefined);
  expect(await page.evaluate(() => window.__cruxCurrent?.length ?? null), 'molecules made the current object').toBe(n);
}, {description: 'read once the sketcher has nothing on its way'});

export const currentMolecule = Then('a molecule should have become the current object {int} time(s), the molecule {string}',
  async (page: Page, n: number, smiles: string) => {
    await expect.poll(() => page.evaluate(() => window.__cruxCurrent?.length ?? null), {message: 'molecules made the current object'}).toBe(n);
    const same = await page.evaluate(async (want) => {
      const last = window.__cruxCurrent![window.__cruxCurrent!.length - 1].value;
      const rdkit = await grok.functions.call('Chem:getRdKitModule');
      const canonical = (v: string) => {
        const mol = rdkit.get_mol(v);
        try {
          return mol.get_smiles();
        } finally {
          mol.delete();
        }
      };
      return {last: canonical(last), want: canonical(want)};
    }, smiles);
    expect(same.last, 'the last molecule made the current object').toBe(same.want);
  }, {description: 'the count, and the last one compared with the molecule named by RDKit'});

// ---------------------------------------------------------------- cells, byte for byte

const cells = new WeakMap<Page, string>();

export const rememberCell = When('user remembers the value of {string} column in row {int}', async (page: Page, column: string, row: number) => {
  cells.set(page, await page.evaluate(([c, r]) => String(grok.shell.t.col(c).get(r - 1) ?? ''), [column, row] as [string, number]));
}, {tier: 'api', description: 'the cell as text, rows counted from 1'});

export const cellAsRemembered = Then('the value of {string} column in row {int} should be, byte for byte, as remembered', async (page: Page, column: string, row: number) => {
  const before = cells.get(page);
  if (before === undefined)
    throw new Error('no cell remembered: "user remembers the value of … column in row …" first');
  await expect.poll(() => page.evaluate(([c, r]) => String(grok.shell.t.col(c).get(r - 1) ?? ''), [column, row] as [string, number]),
    {message: `${column} in row ${row}`}).toBe(before);
}, {description: 'the cell\'s text is the very string remembered'});

// ---------------------------------------------------------------- the sketcher chosen

export const sketcherChosen = Then('the session\'s sketcher and the account\'s choice should be {string}', async (page: Page, name: string) => {
  // the account's choice is what the server holds: the platform sends a setting on a timer of its own, so a page
  // reloaded before it would start from the earlier one
  await expect.poll(() => page.evaluate(async () => [DG.chem.currentSketcherType,
    String(await grok.dapi.userDataStorage.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY) ?? '')]),
  {message: 'the session\'s sketcher type, and the account\'s sketcher setting on the server', timeout: pollMs(30000)}).toEqual([name, name]);
}, {description: 'DG.chem.currentSketcherType, and the account\'s `sketcher/selected` setting as the server holds it (read back until it does)'});

// ---------------------------------------------------------------- the options menu's Recent and Favorites

const COLLECTIONS = ['chem-molecule-recent', 'chem-molecule-favorites'];

export const emptyCollections = Given('the sketcher\'s Recent and Favorites are empty, and come back when the feature ends',
  async (page: Page) => {
    const before = await page.evaluate((keys) => keys.map((k) => localStorage.getItem(k)), COLLECTIONS);
    await page.evaluate((keys) => keys.forEach((k) => localStorage.removeItem(k)), COLLECTIONS);
    atFeatureEnd(page, async () => {
      await page.evaluate(([keys, values]) => keys.forEach((k, i) => values[i] === null ? localStorage.removeItem(k) :
        localStorage.setItem(k, values[i]!)), [COLLECTIONS, before] as [string[], (string | null)[]]);
    });
  }, {tier: 'api', description: 'the molecules the sketcher\'s options menu lists under Recent and Favorites, kept in the browser (DG.chem.Sketcher.RECENT_KEY, FAVORITES_KEY): emptied now, as they were at feature end; the first molecule each gets is written at once'});

export const pickMenuMolecule = When('user picks molecule {int} of the {string} group from the open menu',
  async (page: Page, n: number, group: string) => {
    const item = viewers.menuItems(page, group).filter({visible: true}).first();
    const box = await item.boundingBox();
    if (box === null)
      throw new Error(`the open menu has no "${group}" group`);
    // a Dart group opens on a pointer move into it, as the library's menu picks open one
    const y = box.y + box.height / 2;
    await page.mouse.move(Math.max(0, box.x - 8), y);
    await page.mouse.move(box.x + box.width * 0.75, y);
    await page.mouse.move(box.x + box.width / 2, y);
    const molecules = item.locator(`[name="div-${group}---"]`);
    await expect(molecules.nth(n - 1), `molecule ${n} of the "${group}" group`).toBeVisible();
    const target = await molecules.nth(n - 1).boundingBox();
    await page.mouse.move(target!.x + 4, y);
    await page.mouse.move(target!.x + 4, target!.y + target!.height / 2);
    await molecules.nth(n - 1).click();
  }, {tier: 'ui', description: 'the menu lists a molecule as its drawing, with no label: the nth drawing in the group'});

// ---------------------------------------------------------------- what the clipboard holds

export const clipboardFormat = Then('the clipboard should hold {word} text', async (page: Page, format: string) => {
  if (format !== 'SMILES' && format !== 'MOLBLOCK')
    throw new Error(`a format the clipboard holds is SMILES or MOLBLOCK, not "${format}"`);
  const text = await gestures.readClipboard(page);
  const molblock = text.includes('M  END');
  expect({oneLine: !text.trim().includes('\n'), molblock}, `the clipboard: ${text.slice(0, 80)}`)
    .toEqual(format === 'SMILES' ? {oneLine: true, molblock: false} : {oneLine: false, molblock: true});
}, {description: 'SMILES: one line and no molblock end; MOLBLOCK: a molblock, with its "M  END"'});
