/* Crux Sketch, Chem's own molecule sketcher (src/crux/crux-sketcher.ts), as Chem's features drive it. The sketcher's
   widget ("crux sketcher widget": its "atom N" and "bond N" areas, its readings), its controls by test id and the
   molecule readings are the library's `molecules` tier (bdd.config.json), shared with the other packages whose hosts
   open Crux. Here is what is Chem's own: the drawing a sketcher shows in place of itself (a filter card's, a pane's),
   Ketcher's canvas after a switch, and what the platform does with the sketcher that no gesture shows: the change
   events of the sketcher a step opens next, counted from its creation, and the molecules made the current object. */
import {Page} from '@playwright/test';
import {element, Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, el, expect, gestures, guideFrame, locate, pollMs, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

element('sketcher thumbnail', {selector: '.chem-external-sketcher-canvas',
  description: 'the drawing a sketcher shows in place of itself until clicked (a filter card, a pane\'s scaffold): its Clear button shows while the pointer is over it'});

element('crux aromatic bond tool', {selector: '[data-u2-name="cruxSketch"] crux-sketch [data-testid="toolbar.bond.aromatic"]',
  description: 'the aromatic bond tool of Crux\'s left toolbar: a click on a bond makes it aromatic'});

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

// ---------------------------------------------------------------- the guides about Crux (features/guides/crux/)

/* Crux's own controls by their test ids (crux-sketch's docs/conventions/test-ids.md), beside the molecules tier's
   (crux canvas, crux clear button, crux single bond tool, crux benzene tool, crux nitrogen tool, crux select tool,
   crux undo button, crux flip horizontal button, crux R-group tool and its dialog's buttons): the ones the guides
   about Crux show. Its dialogs and menus are in its open shadow root too. */
const CRUX = '[data-u2-name="cruxSketch"] crux-sketch';
const CRUX_CONTROLS: [string, string, string][] = [
  // drawing
  ['crux double bond tool', 'toolbar.bond.double', 'the double bond tool of Crux\'s left toolbar: a click on a bond makes it double'],
  ['crux dative bond tool', 'toolbar.bond.dative', 'the dative bond tool: a drag draws a dative bond from its donor, where it starts, to its acceptor'],
  ['crux chain tool', 'toolbar.chain', 'the chain tool: a drag from an atom or the empty canvas draws a zig-zag carbon chain, a bond per bond length dragged'],
  ['crux stereo bond tool', 'toolbar.bond.stereo', 'the stereo bonds palette\'s button, the wedge until another member is chosen; a click while its member is the tool opens the palette'],
  ['crux hash bond tool', 'toolbar.bond.hash', 'the hashed wedge, in the open stereo bonds palette'],
  ['crux cyclohexane tool', 'toolbar.ring.cyclohexane', 'the cyclohexane of Crux\'s ring bar: a click on a bond fuses the ring onto it'],
  ['crux cyclopropane tool', 'toolbar.ring.cyclopropane', 'the cyclopropane of Crux\'s ring bar: a click on an atom adds a spiro ring there'],
  ['crux oxygen tool', 'toolbar.element.o', 'the oxygen of Crux\'s element palette: a click on an atom makes it O'],
  ['crux charge plus tool', 'toolbar.charge.plus', 'Charge plus: a click on an atom raises its charge by one'],
  ['crux charge minus tool', 'toolbar.charge.minus', 'Charge minus: a click on an atom lowers its charge by one'],
  ['crux periodic table button', 'toolbar.element.table', 'the periodic table button, last in Crux\'s element palette'],
  ['crux periodic table', 'dialog.periodic-table', 'Crux\'s periodic table: an element\'s button chooses it, Add makes it the tool'],
  ['crux periodic table list button', 'dialog.periodic-table.mode.list', 'List (query mode): the elements chosen next make an atom list'],
  ['crux periodic table not list button', 'dialog.periodic-table.mode.not-list', 'Not list (query mode): the elements chosen next make a NOT list'],
  ['crux periodic table Q button', 'dialog.periodic-table.generic.q', 'Q, a generic atom in the periodic table (query mode): any atom but carbon and hydrogen'],
  ['crux periodic table Add button', 'dialog.periodic-table.add', 'Add: the element, list or generic chosen becomes the tool'],
  ['crux structure library button', 'toolbar.template.library', 'the Structure Library button, last in Crux\'s ring bar'],
  ['crux structure library', 'dialog.structure-library', 'the Structure Library dialog (Ketcher\'s template library): a template\'s card makes it the tool and closes it'],
  ['crux structure library search', 'dialog.structure-library.search', 'the Structure Library\'s search: the groups show only the templates that match'],
  ['crux query bond tool', 'toolbar.bond.query', 'the query bonds palette\'s button (query mode), Any bond until another member is chosen'],
  ['crux attachment point tool', 'toolbar.rgroup.attachment', 'the attachment point tool, in the open R-group palette: a click on an atom opens the Attachment Points dialog'],
  ['crux primary attachment point checkbox', 'dialog.attachment-points.primary', 'Primary attachment point, in the Attachment Points dialog'],
  ['crux Attachment Points OK button', 'dialog.attachment-points.ok', 'OK of the Attachment Points dialog'],
  // the top bar
  ['crux redo button', 'toolbar.redo', 'Redo, beside Undo on Crux\'s top toolbar'],
  ['crux clean up button', 'toolbar.clean', 'Clean Up on Crux\'s top toolbar: tidies bonds and angles, of the selection when there is one'],
  ['crux copy as button', 'toolbar.copy.open', 'Copy As, the opener beside Copy: the formats the drawing is copied in'],
  ['crux SMILES item', 'menu.copy-as.smiles', 'SMILES, in the Copy As menu'],
  ['crux paste button', 'toolbar.paste', 'Paste, on Crux\'s top toolbar'],
  ['crux hydrogens button', 'toolbar.hydrogens', 'Hydrogens, at the end of the top toolbar\'s Structure group: its menu adds or removes explicit hydrogens'],
  ['crux add hydrogens item', 'menu.hydrogens.add', 'Add explicit hydrogens, in the Hydrogens menu'],
  ['crux remove hydrogens item', 'menu.hydrogens.remove', 'Remove explicit hydrogens, in the Hydrogens menu'],
  ['crux settings button', 'toolbar.options', 'the gear, "Settings and help", last on Crux\'s top toolbar'],
  ['crux settings item', 'menu.options.settings', 'Settings…, in the gear\'s menu'],
  ['crux query mode item', 'menu.options.mode', 'Query mode, in the gear\'s menu: switches the sketcher between molecule and query mode'],
  ['crux stereo labels checkbox', 'dialog.settings.stereo-labels', '"Show R, S, E and Z labels", in Crux\'s Settings'],
  ['crux settings Apply button', 'dialog.settings.apply', 'Apply of Crux\'s Settings'],
  // the canvas and its menus
  ['crux rotate handle', 'overlay.rotate-handle', 'the handle above a selection of two atoms or more: a drag of it turns the selection'],
  ['crux formula', 'panel.info.formula', 'the formula in Crux\'s formula and mass readout, at the canvas\'s foot'],
  ['crux context menu', 'menu.context', 'Crux\'s own context menu, of an atom, a bond, the selection or the canvas'],
  ['crux Atom Properties item', 'menu.context.atom-properties', 'Atom Properties…, in an atom\'s context menu'],
  ['crux hash item', 'menu.context.bond.hash', 'the hashed wedge among the bond types of a bond\'s context menu'],
  ['crux wedge item', 'menu.context.bond.wedge', 'the wedge among the bond types of a bond\'s context menu'],
  ['crux atropisomer P item', 'menu.context.atrop.p', 'P, in an atropisomer axis bond\'s context menu'],
  ['crux atropisomer M item', 'menu.context.atrop.m', 'M, in an atropisomer axis bond\'s context menu'],
  ['crux E item', 'menu.context.ez.e', 'E, in a stereo double bond\'s context menu'],
  ['crux Z item', 'menu.context.ez.z', 'Z, in a stereo double bond\'s context menu'],
  ['crux enhanced stereo item', 'menu.context.enhanced-stereo', 'Enhanced Stereochemistry…, in the context menu of a stereocentre or the selection'],
  ['crux topology item', 'menu.context.topology', 'Topology, in a bond\'s context menu (query mode): opens its submenu'],
  ['crux ring topology item', 'menu.context.topology.ring', 'Ring, in the Topology submenu'],
  ['crux chain topology item', 'menu.context.topology.chain', 'Chain, in the Topology submenu'],
  ['crux expand abbreviation item', 'menu.context.expand-abbreviation', 'Expand Abbreviation, in a contracted abbreviation\'s context menu'],
  ['crux contract abbreviation item', 'menu.context.contract-abbreviation', 'Contract Abbreviation, in the context menu of an expanded abbreviation\'s atom'],
  // dialogs
  ['crux Atom Properties', 'dialog.atom-properties', 'Crux\'s Atom Properties dialog'],
  ['crux charge field', 'dialog.atom-properties.charge', 'Charge, in Atom Properties'],
  ['crux isotope field', 'dialog.atom-properties.isotope', 'Isotope, in Atom Properties'],
  ['crux radical list', 'dialog.atom-properties.radical', 'Radical, in Atom Properties: none, a monovalent or a divalent radical'],
  ['crux Atom Properties Apply button', 'dialog.atom-properties.apply', 'Apply of Atom Properties'],
  ['crux enhanced stereo dialog', 'dialog.enhanced-stereo', 'Crux\'s Enhanced Stereochemistry dialog'],
  ['crux AND button', 'dialog.enhanced-stereo.type.and', 'AND, in Enhanced Stereochemistry'],
  ['crux OR button', 'dialog.enhanced-stereo.type.or', 'OR, in Enhanced Stereochemistry'],
  ['crux enhanced stereo Apply button', 'dialog.enhanced-stereo.apply', 'Apply of Enhanced Stereochemistry'],
];
for (const [name, id, description] of CRUX_CONTROLS)
  element(name, {selector: `${CRUX} [data-testid="${id}"]`, description});

/** The sketcher a guide about Crux opens on, sized for the video, in the middle of the page; its zoom, in Zoom in
 * steps from 100%, so the drawing reads in a video. */
const GUIDE_DIALOG = {width: 1040, height: 860, zoomSteps: 5};

/** Crux's settings as a user leaves them, kept in the page (crux-sketch's `persist-settings`). */
const CRUX_SETTINGS = 'crux-sketch:settings';

/** In the page: Crux's status values (its `getWidgetStatus()` readings). */
function cruxValues(page: Page): Promise<any> {
  return viewers.onViewer(page, el('crux sketcher widget'),
    (e) => (window as any).__bdd.viewerOf(e).getWidgetStatus()?.values ?? {});
}

async function openCrux(page: Page, molecule: string, labels: boolean): Promise<void> {
  const {width, height, zoomSteps} = GUIDE_DIALOG;
  // set-up a person never does: the video starts with the sketcher open
  silent(page);
  const kept = await page.evaluate((k) => localStorage.getItem(k), CRUX_SETTINGS);
  atFeatureEnd(page, async () => {
    await page.evaluate(([k, v]) => v === null ? localStorage.removeItem(k) : localStorage.setItem(k, v), [CRUX_SETTINGS, kept]);
  });
  await page.evaluate(([m, w, h]) => {
    const sketcher = new DG.chem.Sketcher();
    if (m)
      sketcher.setMolecule(m);
    (window as any).ui.dialog().add(sketcher).onOK(() => {})
      .show({resizable: true, width: w, height: h, x: Math.round((innerWidth - w) / 2), y: Math.round((innerHeight - h) / 2)});
  }, [molecule, width, height] as [string, number, number]);
  await expect.poll(async () => {
    const v = await cruxValues(page).catch(() => ({}));
    return v.ready === true && v.pending === false && (molecule === '' || Number(v.atoms) > 0);
  }, {message: 'Crux ready in the sketcher dialog, showing the molecule'}).toBe(true);
  if (labels)
    await page.evaluate((root) => { (document.querySelector(root) as any).stereoLabels = 'shown'; }, CRUX);
  // Zoom in (F7) on the canvas, as a person makes a drawing larger
  await page.evaluate((root) => (document.querySelector(root)?.shadowRoot?.querySelector('[data-testid="canvas"]') as HTMLElement | null)
    ?.focus(), CRUX);
  for (let i = 0; i < zoomSteps; i++)
    await page.keyboard.press('F7');
  await viewers.settle(page, el('crux sketcher widget'));
  const dialog = (await locate(page, el('sketcher dialog'))).filter({visible: true}).first();
  const box = await dialog.boundingBox();
  // the menus and dialogs Crux opens over the canvas can reach past the dialog's edges
  if (box)
    guideFrame(page, box, 40);
}

export const cruxOpenOn = Given('the Crux sketcher is open on {string}', (page: Page, molecule: string) => openCrux(page, molecule, false),
  {tier: 'api', description: 'a sketcher dialog with Crux (pinned before) showing the molecule, sized for a video, which a guide shows alone; Crux\'s settings come back at feature end; silent in a guide'});

export const cruxOpenEmpty = Given('the Crux sketcher is open', (page: Page) => openCrux(page, '', false),
  {tier: 'api', description: 'the same, empty'});

export const cruxOpenLabelled = Given('the Crux sketcher is open on {string}, showing R, S, E and Z labels', (page: Page, molecule: string) =>
  openCrux(page, molecule, true),
{tier: 'api', description: 'the same, with Crux\'s CIP labels shown (its stereoLabels), as Settings\' "Show R, S, E and Z labels" shows them'});

export const cruxOpenOnMolfile = Given('the Crux sketcher is open on this molfile, showing R, S, E and Z labels:',
  (page: Page, molfile: string) => openCrux(page, molfile, true),
  {tier: 'api', description: 'the same, on the molfile the doc string holds (a drawing with its wedges where they were drawn)'});

const STEREO_DRAWN: {[stereo: string]: string} = {wedge: 'up', hash: 'down', plain: 'none'};

export const cruxBondDrawn = Then('bond {int} of Crux should be drawn as a {word}', async (page: Page, bond: number, kind: string) => {
  const want = STEREO_DRAWN[kind];
  if (want === undefined)
    throw new Error(`a bond is drawn as a wedge, a hash or plain, not "${kind}"`);
  await expect(page.locator(`${CRUX} [data-testid="bond.${bond}"]`), `bond ${bond}'s stereo as drawn`).toHaveAttribute('data-stereo', want);
}, {description: 'the bond\'s data-stereo as Crux draws it: wedge (up), hash (down) or plain (none)'});

export const cruxKey = When('user presses the {string} key over the {string} area of {widget}',
  async (page: Page, key: string, area: string, target: ElementRef) => {
    const c = viewers.centerOf(await viewers.hitArea(page, target, area, true));
    await page.mouse.move(c.x - 3, c.y - 3);
    await page.mouse.move(c.x, c.y);
    await viewers.settle(page, target);
    // a key goes to the sketcher holding the focus: a person has clicked in it before
    await page.evaluate((root) => {
      const sketch = document.querySelector(root) as HTMLElement | null;
      if (sketch && !sketch.contains(document.activeElement) && sketch.shadowRoot?.activeElement == null)
        (sketch.shadowRoot?.querySelector('[data-testid="canvas"]') as HTMLElement | null)?.focus();
    }, CRUX);
    // as typed: Crux's keymap is case-sensitive (o is oxygen, Shift+O a methoxy)
    await page.keyboard.press(key);
  }, {tier: 'ui', description: 'the pointer over the area, then the key as written, case and all (Crux\'s hotkeys act on the hovered atom or bond)'});

export const cruxMenu = When('user opens the Crux context menu on the {string} area', async (page: Page, area: string) => {
  const c = viewers.centerOf(await viewers.hitArea(page, el('crux sketcher widget'), area, true));
  await page.mouse.click(c.x, c.y, {button: 'right'});
  await expect(page.locator(`${CRUX} [data-testid="menu.context"]`), 'Crux\'s context menu').toBeVisible();
}, {tier: 'ui', description: 'a right-click on an atom or a bond of Crux (an area of crux sketcher widget): Crux\'s own menu, in its shadow root'});

export const cruxSpot = When('user clicks on Crux canvas {int}% across and {int}% down', async (page: Page, across: number, down: number) => {
  const b = await viewers.hitArea(page, el('crux sketcher widget'), 'canvas', true);
  await page.mouse.click(b.x + b.width * across / 100, b.y + b.height * down / 100);
}, {tier: 'ui', description: 'a click on a point of Crux\'s canvas, by its place across and down the canvas (an empty spot)'});

export const cruxSmarts = Then('the Crux sketcher should hold the query {string}', async (page: Page, smarts: string) => {
  await viewers.settle(page, el('crux sketcher widget')).catch(() => undefined);
  await expect.poll(() => viewers.onViewer(page, el('crux sketcher widget'),
    (e) => (window as any).DG.Widget.find(e).getSmarts()), {message: 'the SMARTS of Crux\'s drawing'}).toBe(smarts);
}, {description: 'the SMARTS Crux writes for its drawing (getSmarts), as written'});
