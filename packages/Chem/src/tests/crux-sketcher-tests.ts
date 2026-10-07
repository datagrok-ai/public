import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

const SMILES = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const STEREO_SMILES = 'N[C@@H](C)C(=O)O';
const SMARTS = '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
// The same ring with an atom list: a query feature, which no SMILES holds
const SMARTS_LIST = '[#6,#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
const SMARTS_TARGET = 'Cc1ccc(O)cc1';
/** Long enough for the host's handlers and for any further event a change could still cause. */
const SETTLE_MS = 300;

interface OpenSketcher {
  host: DG.chem.Sketcher;
  dialog: DG.Dialog;
  /** How many times the implementation's onChanged fired since it was created. */
  changes: () => number;
}

interface CruxElement extends HTMLElement {
  readonly isEmpty: boolean;
  readonly mode: string;
  readonly molfile: string;
  readonly smarts: string;
  /** Where Crux draws each atom and bond (API-056), in px from the element's corner. */
  readonly positions: {atoms: ({x: number, y: number} | null)[], bonds: ({x: number, y: number} | null)[]};
}

/** The getters of an implementation read at one moment, and what reading one threw. */
interface Reads {
  smiles?: string;
  molFile?: string;
  molV3000?: string;
  threw?: string;
}

/** Molecules a SMARTS is matched over, so that two queries are compared by what they find. */
const PROBES = ['c1ccccc1', 'c1ccncc1', 'Cc1ccccc1', 'c1ccc2ccccc2c1', 'CC(=O)Nc1ccc(O)cc1', 'CC(=O)OC', 'CCN',
  'OCC(=O)N', 'CC', 'CCC(=O)O', 'C1CCCCC1'];

let rdkit: any;

function canonical(molecule: string): string {
  const mol = rdkit.get_mol(molecule);
  if (!mol)
    throw new Error(`RDKit cannot read: ${molecule}`);
  try {
    return mol.get_smiles();
  } finally {
    mol.delete();
  }
}

function molblocks(smiles: string): {v2000: string, v3000: string} {
  const mol = rdkit.get_mol(smiles);
  try {
    return {v2000: mol.get_molblock(), v3000: mol.get_v3Kmolblock()};
  } finally {
    mol.delete();
  }
}

function rdkitSmarts(smiles: string): string {
  const mol = rdkit.get_mol(smiles);
  try {
    return mol.get_smarts();
  } finally {
    mol.delete();
  }
}

/** The atom sets of `target` that `query` (SMARTS or a molblock) matches, as RDKit finds them. */
function matchSets(target: string, query: string): string {
  const mol = rdkit.get_mol(target);
  const qmol = rdkit.get_qmol(query);
  try {
    if (!qmol)
      throw new Error(`RDKit cannot read the query: ${query}`);
    const found = JSON.parse(mol.get_substruct_matches(qmol) || '[]');
    const sets = (Array.isArray(found) ? found : []).map((m: {atoms: number[]}) => [...m.atoms].sort((x, y) => x - y));
    return JSON.stringify(sets.sort());
  } finally {
    qmol?.delete();
    mol.delete();
  }
}

async function openCrux(validation?: (s: string) => string | null, substructureFilter = false): Promise<OpenSketcher> {
  return openWith({validation, substructureFilter});
}

interface OpenOptions {
  validation?: (s: string) => string | null;
  substructureFilter?: boolean;
  /** What the host shows as the sketcher opens: SMILES, or SMARTS. */
  smiles?: string;
  smarts?: string;
  /** The dialog's OK handler, as a host application gives it. */
  onOK?: () => void;
}

/** Crux in DG.chem.Sketcher in a dialog, as Datagrok's own code opens a sketcher: its change events counted from its
 * creation, so a sketcher opened showing a molecule has fired once when this resolves. */
async function openWith(o: OpenOptions): Promise<OpenSketcher> {
  grok.chem.currentSketcherType = 'Crux';
  const host = new DG.chem.Sketcher(undefined, o.validation);
  // As Chem's substructure filter does, before the sketcher is created
  if (o.substructureFilter)
    host.isSubstructureFilter = true;
  if (o.smiles !== undefined)
    host.setSmiles(o.smiles);
  if (o.smarts !== undefined)
    host.setSmarts(o.smarts);
  let count = 0;
  let created: DG.chem.SketcherBase | null = null;
  // the host puts the sketcher's root in its box before it awaits its init: its events are counted from then
  const observer = new MutationObserver(() => {
    if (host.sketcher === null || host.sketcher === created)
      return;
    created = host.sketcher;
    created.onChanged.subscribe(() => count++);
  });
  observer.observe(host.host, {childList: true});
  const dialog = ui.dialog().add(host);
  if (o.onOK)
    dialog.onOK(o.onOK);
  dialog.show();
  try {
    await host.sketcherReady();
  } finally {
    observer.disconnect();
  }
  expect(created !== null && created === host.sketcher, true, 'the sketcher was not seen as it was created');
  expect(host.sketcher!.root.classList.contains('crux-sketcher'), true, 'the sketcher is not Crux');
  await delay(SETTLE_MS);
  return {host, dialog, changes: () => count};
}

function cruxElement(s: OpenSketcher): CruxElement | null {
  return s.host.sketcher!.root.querySelector<CruxElement>('crux-sketch');
}

/** Presses a toolbar button of the sketcher, as the user does. */
function press(s: OpenSketcher, tool: string): void {
  const button = cruxElement(s)!.shadowRoot!.querySelector<HTMLButtonElement>(`[data-testid="toolbar.${tool}"]`);
  if (!button)
    throw new Error(`no toolbar button ${tool}`);
  button.click();
}

/** A user's press and release on Crux's canvas at a point of the element (its `positions`' frame), as the pointer's
 * events reach the canvas. */
function pressAt(s: OpenSketcher, p: {x: number, y: number}): void {
  const el = cruxElement(s)!;
  const canvas = el.shadowRoot!.querySelector('[data-testid="canvas"]')!;
  const r = el.getBoundingClientRect();
  const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true,
    clientX: r.left + p.x, clientY: r.top + p.y, pointerId: 1, pointerType: 'mouse', isPrimary: true, button: 0,
    buttons});
  canvas.dispatchEvent(new PointerEvent('pointerdown', init(1)));
  canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
}

/** Where Crux draws atom `i` (its `positions`). */
function atomAt(s: OpenSketcher, i: number): {x: number, y: number} {
  const p = cruxElement(s)!.positions.atoms[i];
  if (!p)
    throw new Error(`Crux draws no atom ${i}`);
  return p;
}

/** Crux's canvas's box, in px from the element's corner. */
function canvasBox(s: OpenSketcher): {x: number, y: number, width: number, height: number} {
  const el = cruxElement(s)!;
  const c = el.shadowRoot!.querySelector('[data-testid="canvas"]')!.getBoundingClientRect();
  const r = el.getBoundingClientRect();
  return {x: c.left - r.left, y: c.top - r.top, width: c.width, height: c.height};
}

/** The middle of Crux's canvas. */
function canvasMiddle(s: OpenSketcher): {x: number, y: number} {
  const c = canvasBox(s);
  return {x: c.x + c.width / 2, y: c.y + c.height / 2};
}

function readAll(impl: DG.chem.SketcherBase): Reads {
  const out: Reads = {};
  for (const k of ['smiles', 'molFile', 'molV3000'] as const) {
    try {
      out[k] = impl[k];
    } catch (e) {
      out.threw = `${k}: ${e}`;
    }
  }
  return out;
}

/** The three getters read `smiles`'s molecule, the molblocks V2000 and V3000. */
function expectReads(got: Reads | null, smiles: string, where: string): void {
  expect(got !== null, true, `${where}: nothing read`);
  expect(got!.threw ?? null, null, `${where}: a getter threw`);
  const want = canonical(smiles);
  expect(canonical(got!.smiles!), want, `${where}: smiles`);
  expect(canonical(got!.molFile!), want, `${where}: molFile`);
  expect(got!.molFile!.includes('V2000'), true, `${where}: molFile is V2000`);
  expect(canonical(got!.molV3000!), want, `${where}: molV3000`);
  expect(got!.molV3000!.includes('V3000'), true, `${where}: molV3000 is V3000`);
}

function atomCount(smiles: string): number {
  const mol = rdkit.get_mol(smiles);
  try {
    return JSON.parse(mol.get_json()).molecules[0].atoms.length;
  } finally {
    mol.delete();
  }
}

/** What Crux itself holds (its element, not the adapter's getters) is `smiles`'s molecule, and it draws its atoms. */
function expectDrawing(el: CruxElement | null, smiles: string): void {
  expect(el !== null, true, 'no <crux-sketch>');
  expect(el!.isEmpty, false, `Crux's canvas, for ${smiles}`);
  expect(canonical(el!.molfile), canonical(smiles), 'the molecule Crux holds');
  expect(el!.positions.atoms.length, atomCount(smiles), 'the atoms Crux draws');
}

function expectEmptyCanvas(el: CruxElement | null): void {
  expect(el !== null, true, 'no <crux-sketch>');
  expect(el!.isEmpty, true, 'Crux\'s canvas is empty');
  expect(el!.positions.atoms.length, 0, 'the atoms Crux draws');
}

/** The SMARTS `a` finds what `b` finds, over the probes. */
function sameQuery(a: string, b: string): boolean {
  return PROBES.every((p) => matchSets(p, a) === matchSets(p, b));
}

/** Opens the host's options menu (≡) by its icon, as the user does, and finds the item named `name` in it. */
async function optionsMenuItem(host: DG.chem.Sketcher, name: string): Promise<HTMLElement> {
  host.root.querySelector<HTMLElement>('.grok-sketcher-input .d4-input-options')!.click();
  let item: HTMLElement | null = null;
  await awaitCheck(() => (item = document.querySelector<HTMLElement>(`.d4-menu-popup [d4-name="${name}"]`)) !== null,
    `no "${name}" in the options menu`, 5000);
  return item!;
}

/** Runs a change and checks that it fires the implementation's onChanged exactly once. */
async function expectOneChange(s: OpenSketcher, change: () => void): Promise<void> {
  const before = s.changes();
  change();
  await delay(SETTLE_MS);
  expect(s.changes() - before, 1, 'onChanged events');
}

category('Crux sketcher', () => {
  let previousSketcherType: string;

  before(async () => {
    previousSketcherType = grok.chem.currentSketcherType;
    rdkit = await grok.functions.call('Chem:getRdKitModule');
  });

  after(async () => {
    grok.chem.currentSketcherType = previousSketcherType;
  });

  test('registered', async () => {
    const funcs = DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}});
    expect(funcs.some((f) => f.friendlyName === 'Crux'), true, 'Crux is not a molecule sketcher');
    // HOST-001, HOST-003: registered by Chem alone; a choice of Chem's Sketcher setting, never its default
    expect(funcs.filter((f) => f.friendlyName === 'Crux').map((f) => f.package?.name).join(', '), 'Chem',
      'the packages that register Crux');
    const sketcher = (await (await fetch(`${_package.webRoot}package.json`)).json())
      .properties.find((p: {name: string}) => p.name === 'Sketcher');
    expect(sketcher.choices.includes('Crux'), true, `Chem's Sketcher setting offers ${sketcher.choices.join(', ')}`);
    expect(sketcher.defaultValue, 'Ketcher', 'the default of Chem\'s Sketcher setting');
    const s = await openCrux();
    try {
      expect(s.changes(), 0, 'onChanged events on open');
      expect(cruxElement(s) !== null, true, 'no <crux-sketch> in the sketcher');
      expect(cruxElement(s)!.mode, 'molecule', 'the mode in a host that asks for a molecule');
    } finally {
      s.dialog.close();
    }
  });

  test('query mode in the substructure filter', async () => {
    const s = await openCrux(undefined, true);
    try {
      expect(cruxElement(s)!.mode, 'query', 'the mode in the substructure filter');
      expect(s.changes(), 0, 'onChanged events on open');
    } finally {
      s.dialog.close();
    }
  });

  test('size', async () => {
    const s = await openCrux();
    try {
      const sketcher = s.host.sketcher!;
      expect(`${sketcher.width} x ${sketcher.height}`, '0 x 0', 'the minimum size Crux reports');
      for (const e of [s.host.host, sketcher.root, cruxElement(s)!]) {
        const style = getComputedStyle(e);
        expect(['0px', 'auto'].includes(style.minWidth) && ['0px', 'auto'].includes(style.minHeight), true,
          `min size of ${e.className || e.tagName}: ${style.minWidth} x ${style.minHeight}`);
      }
      const size = (e: Element) => {
        const r = e.getBoundingClientRect();
        return `${Math.round(r.width)} x ${Math.round(r.height)}`;
      };
      // a dialog opens at the sketcher's own size; Crux fills the box when the user shrinks or grows the dialog
      expect(size(cruxElement(s)!), '500 x 400', 'the size a dialog opens at');
      for (const [w, h] of [[300, 330], [200, 260], [900, 760]]) {
        s.dialog.root.style.width = `${w}px`;
        s.dialog.root.style.height = `${h}px`;
        await delay(SETTLE_MS);
        const box = s.host.host.getBoundingClientRect();
        expect(box.width < w && box.height < h, true, `the host box ${size(s.host.host)} in a ${w} x ${h} dialog`);
        expect(size(cruxElement(s)!), size(s.host.host), `Crux in the host box of a ${w} x ${h} dialog`);
      }
    } finally {
      s.dialog.close();
    }
  });

  test('smiles', async () => {
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      expect(s.host.getSmiles(), SMILES, 'explicitMol');
      expect(canonical(s.host.getMolFile()), canonical(SMILES), 'V2000 written by Crux');
      expect(canonical(s.host.sketcher!.molV3000), canonical(SMILES), 'V3000 written by Crux');
      expect(matchSets(SMILES, (await s.host.getSmarts())!), matchSets(SMILES, rdkitSmarts(SMILES)),
        'SMARTS written by Crux');
      await expectOneChange(s, () => s.host.setSmiles(STEREO_SMILES));
      expect(canonical(s.host.getMolFile()), canonical(STEREO_SMILES), 'stereo in the V2000 written by Crux');
    } finally {
      s.dialog.close();
    }
  });

  test('molV2000', async () => {
    const {v2000} = molblocks(SMILES);
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setMolFile(v2000));
      expect(s.host.getMolFile(), v2000, 'explicitMol');
      expect(canonical(s.host.getSmiles()), canonical(SMILES), 'SMILES written by Crux');
      expect(s.host.sketcher!.molV3000.includes('V3000'), true, 'V3000 written by Crux');
      expect(canonical(s.host.sketcher!.molV3000), canonical(SMILES), 'V3000 written by Crux');
    } finally {
      s.dialog.close();
    }
  });

  test('molV3000', async () => {
    const {v3000} = molblocks(SMILES);
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setMolFile(v3000));
      expect(s.host.getMolFile(), v3000, 'explicitMol');
      expect(canonical(s.host.getSmiles()), canonical(SMILES), 'SMILES written by Crux');
      expect(s.host.sketcher!.molFile.includes('V2000'), true, 'V2000 written by Crux');
      expect(canonical(s.host.sketcher!.molFile), canonical(SMILES), 'V2000 written by Crux');
    } finally {
      s.dialog.close();
    }
  });

  test('smarts', async () => {
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setSmarts(SMARTS));
      expect(await s.host.getSmarts(), SMARTS, 'explicitMol');
      // the substructure filter searches with the SMARTS that RDKit makes of the V2000 the sketcher writes
      const smarts = DG.chem.convert(s.host.getMolFile(), DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
      expect(matchSets(SMARTS_TARGET, smarts), matchSets(SMARTS_TARGET, SMARTS), 'query V2000 written by Crux');
      // plain atoms and aromatic bonds: Crux writes the SMILES RDKit writes of the same molecule
      expect(canonical(s.host.sketcher!.smiles), canonical('c1ccccc1'), 'SMILES of a SMARTS without query features');
      await expectOneChange(s, () => s.host.setSmarts(SMARTS_LIST));
      expect(s.host.sketcher!.smiles, '', 'SMILES of a query');
      const listSmarts = DG.chem.convert(s.host.getMolFile(), DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
      expect(matchSets('Cc1ccncc1', listSmarts), matchSets('Cc1ccncc1', SMARTS_LIST), 'atom list V2000 written by Crux');
    } finally {
      s.dialog.close();
    }
  });

  test('value set while loading', async () => {
    grok.chem.currentSketcherType = 'Crux';
    const host = new DG.chem.Sketcher();
    let count = 0;
    let loading: boolean | null = null;
    // the host puts the new sketcher's root in its box and then awaits its init: a value set now waits for Crux
    const observer = new MutationObserver(() => {
      if (host.sketcher === null || loading !== null)
        return;
      host.sketcher.onChanged.subscribe(() => count++);
      loading = !host.sketcher.isInitialized;
      host.setSmarts(SMARTS);
    });
    observer.observe(host.host, {childList: true});
    const dialog = ui.dialog().add(host).show();
    try {
      await host.sketcherReady();
      await delay(SETTLE_MS);
      expect(loading, true, 'the value was set while Crux was loading');
      expect(count, 1, 'onChanged events');
      expect(await host.getSmarts(), SMARTS, 'the value set while loading');
      expect(cruxElement({host, dialog, changes: () => count})!.isEmpty, false, 'canvas');
    } finally {
      observer.disconnect();
      dialog.close();
    }
  });

  test('explicitMol until a user edit', async () => {
    const s = await openCrux();
    try {
      const smiles = 'C1=CC=CC=C1';
      await expectOneChange(s, () => s.host.setSmiles(smiles));
      expect(s.host.sketcher!.explicitMol?.value, smiles, 'explicitMol after the write');
      expect(s.host.getSmiles(), smiles, 'the caller\'s SMILES');
      await expectOneChange(s, () => press(s, 'clear'));
      expect(s.host.sketcher!.explicitMol ?? null, null, 'explicitMol after a user edit');
      expect(s.host.getSmiles(), '', 'SMILES after the user cleared the canvas');
      expect(s.host.isEmpty(), true, 'isEmpty after the user cleared the canvas');
    } finally {
      s.dialog.close();
    }
  });

  test('empty values clear', async () => {
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      await expectOneChange(s, () => s.host.setMolecule(''));
      expect(cruxElement(s)!.isEmpty, true, 'canvas after \'\'');
      expect(s.host.isEmpty(), true, 'isEmpty after \'\'');
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      await expectOneChange(s, () => s.host.setMolFile(DG.WHITE_MOLBLOCK));
      expect(cruxElement(s)!.isEmpty, true, 'canvas after WHITE_MOLBLOCK');
      expect(s.host.isEmpty(), true, 'isEmpty after WHITE_MOLBLOCK');
      s.host.sketcher!.explicitMol = null;
      expect(DG.chem.Sketcher.isEmptyMolfile(s.host.sketcher!.molFile), true, 'the empty canvas\'s V2000');
    } finally {
      s.dialog.close();
    }
  });

  test('unreadable value', async () => {
    const validate = (molecule: string) => DG.Func.find({package: 'Chem', name: 'validateMolecule'})[0]
      .prepare({s: molecule}).callSync().getOutputParamValue();
    const s = await openCrux(validate);
    try {
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      const garbage = 'C1CC(((';
      const changes = s.changes();
      s.host.setSmiles(garbage);
      await delay(SETTLE_MS);
      // as with the other sketchers, the drawing stays and the host keeps the value shown as malformed
      expect(s.changes() - changes, 0, 'onChanged events of an unreadable value');
      expect(canonical(s.host.sketcher!.molFile), canonical(SMILES), 'the drawing after an unreadable value');
      expect(s.host.getSmiles(), garbage, 'explicitMol of an unreadable value');
      expect(s.dialog.root.querySelector('.chem-invalid-molecule-warning')?.textContent, 'Malformed molecule',
        'the host\'s warning');
    } finally {
      s.dialog.close();
    }
  });

  test('two sketchers', async () => {
    const a = await openCrux();
    const b = await openCrux();
    try {
      await expectOneChange(a, () => a.host.setSmiles(SMILES));
      await expectOneChange(b, () => b.host.setSmiles(STEREO_SMILES));
      expect(a.changes(), 1, 'onChanged events of the first sketcher');
      expect(canonical(a.host.sketcher!.molFile), canonical(SMILES), 'the first sketcher\'s molecule');
      expect(canonical(b.host.sketcher!.molFile), canonical(STEREO_SMILES), 'the second sketcher\'s molecule');
      await expectOneChange(a, () => press(a, 'clear'));
      expect(cruxElement(a)!.isEmpty, true, 'the first sketcher after a clear');
      expect(cruxElement(b)!.isEmpty, false, 'the second sketcher after the first was cleared');
      expect(b.changes(), 1, 'onChanged events of the second sketcher');
    } finally {
      a.dialog.close();
      b.dialog.close();
    }
  });

  test('detach', async () => {
    const live = await openCrux();
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      await expectOneChange(s, () => press(s, 'clear'));
      await expectOneChange(s, () => s.host.setSmiles(STEREO_SMILES));
      s.dialog.close();
      await awaitCheck(() => s.host.sketcher!.isDetached, 'Crux has not been detached', 3000);
      expect(cruxElement(s), null, '<crux-sketch> left in the sketcher after detach');
      expect(canonical(s.host.getSmiles()), canonical(STEREO_SMILES), 'SMILES after detach');
      expect(canonical(s.host.sketcher!.molFile), canonical(STEREO_SMILES), 'V2000 after detach');
      expect(canonical(s.host.sketcher!.molV3000), canonical(STEREO_SMILES), 'V3000 after detach');
      s.host.setMolecule('');
      expect(s.host.isEmpty(), true, 'isEmpty after a clear while detached');
      await expectOneChange(live, () => live.host.setSmiles(SMILES));
      expect(canonical(live.host.sketcher!.molFile), canonical(SMILES), 'the other sketcher after detach');
    } finally {
      live.dialog.close();
    }
  });

  test('switch to OpenChemLib and back', async () => {
    const s = await openCrux();
    try {
      await expectOneChange(s, () => s.host.setSmiles(SMILES));
      const crux = s.host.sketcher;
      grok.chem.currentSketcherType = 'OpenChemLib';
      s.host.sketcherType = 'OpenChemLib';
      await s.host.sketcherReady();
      expect(s.host.sketcher !== crux, true, 'OpenChemLib replaced Crux');
      await awaitCheck(() => canonical(s.host.getMolFile()) === canonical(SMILES), 'the molecule in OpenChemLib', 5000);
      const ocl = s.host.sketcher;
      grok.chem.currentSketcherType = 'Crux';
      s.host.sketcherType = 'Crux';
      await s.host.sketcherReady();
      expect(s.host.sketcher !== ocl, true, 'Crux replaced OpenChemLib');
      expect(s.host.sketcher!.root.classList.contains('crux-sketcher'), true, 'the sketcher is not Crux');
      await awaitCheck(() => cruxElement(s)?.isEmpty === false, 'the molecule did not come back to Crux', 5000);
      expect(canonical(s.host.sketcher!.molFile), canonical(SMILES), 'the molecule back in Crux');
    } finally {
      s.dialog.close();
    }
  });
  // HOST-004, its SMARTS row: each switch hands the next sketcher the caller's string, as the options menu (≡) makes a
  // switch (the session's sketcher type, then the host's); the menu itself is crux-choice.feature's
  test('a switch keeps a SMARTS: OpenChemLib, then Crux, Ketcher and Crux', async () => {
    const query = '[#6]1:[#6]:[#6,#7]:[#6]:[#6]:[#6]:1';
    const errors: string[] = [];
    const consoleError = console.error;
    console.error = (...args: unknown[]) => {
      errors.push(args.map(String).join(' '));
      consoleError(...args);
    };
    grok.chem.currentSketcherType = 'OpenChemLib';
    const host = new DG.chem.Sketcher();
    host.setMolecule(query, true);
    const dialog = ui.dialog().add(host).show();
    try {
      await host.sketcherReady();
      for (const type of ['Crux', 'Ketcher', 'Crux']) {
        grok.chem.currentSketcherType = type;
        host.sketcherType = type;
        const sketcher = await host.sketcherReady();
        expect(sketcher.root.classList.contains('crux-sketcher'), type === 'Crux',
          `the sketcher after the switch to ${type}`);
        const smarts = await host.getSmarts();
        expect(sameQuery(smarts!, query), true, `${type}: the host's SMARTS ${smarts} finds what ${query} finds`);
        if (type === 'Crux') {
          const el = sketcher.root.querySelector<CruxElement>('crux-sketch')!;
          expect(!el.isEmpty && sameQuery(el.smarts, query), true, `Crux's own SMARTS ${el.smarts}`);
        }
      }
      expect(errors.filter((e) => /sketch|crux|ketcher|openchemlib|ocl/i.test(e)).join('; '), '',
        'sketcher errors logged');
    } finally {
      console.error = consoleError;
      dialog.close();
    }
  });

  // HOST-002 (spike datagrok-platform, R7b): the options menu (≡) is js-api's. The user's pick is the account's sketcher,
  // which the platform otherwise sends to the server on its 10 s timer: a page reloaded before it lost the pick.
  test('the options menu sends the sketcher picked to the server at once', async () => {
    const saved = grok.userSettings.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY);
    const pick = saved === 'Crux' ? 'OpenChemLib' : 'Crux';
    grok.chem.currentSketcherType = pick === 'Crux' ? 'OpenChemLib' : 'Crux';
    const host = new DG.chem.Sketcher();
    const dialog = ui.dialog().add(host).show();
    const sent: number[] = [];
    const requests = new PerformanceObserver((list) => {
      for (const e of list.getEntries()) {
        if (/\/user_settings_storage\/sketcher\b/.test(e.name))
          sent.push(e.startTime);
      }
    });
    try {
      await host.sketcherReady();
      const item = await optionsMenuItem(host, pick);
      requests.observe({type: 'resource'});
      const picked = performance.now();
      item.click();
      await awaitCheck(() => sent.length > 0, 'the pick never reached the server', 15000);
      const after = Math.round(sent[0] - picked);
      expect(after < 500, true, `the pick was sent ${after} ms after it was made`);
      expect(await grok.dapi.userDataStorage.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY, true), pick,
        'the account\'s sketcher on the server');
      await host.sketcherReady();
    } finally {
      requests.disconnect();
      dialog.close();
      // the account's own choice back on the server
      if (saved === undefined || saved === null)
        grok.userSettings.delete(DG.chem.STORAGE_NAME, DG.chem.KEY);
      else
        grok.userSettings.add(DG.chem.STORAGE_NAME, DG.chem.KEY, saved);
      await grok.userSettings.flush();
    }
  });

  // HOST-002 (R7b): a host keeps its sketcher when the session's changes (Bio's monomer manager keeps its editor's for the
  // page's life); its menu checked the session's sketcher and ignored a pick of it, so that host could not be switched
  // to it.
  test('the options menu checks and switches the sketcher its host shows, not the session\'s', async () => {
    grok.chem.currentSketcherType = 'OpenChemLib';
    const host = new DG.chem.Sketcher();
    const dialog = ui.dialog().add(host).show();
    try {
      const shown = await host.sketcherReady();
      expect(shown.root.classList.contains('crux-sketcher'), false, 'the host opens OpenChemLib');
      // another host's pick, or a run's pin: the session's sketcher is Crux now
      grok.chem.currentSketcherType = 'Crux';
      const item = await optionsMenuItem(host, 'Crux');
      const checked = Array.from(document.querySelectorAll('.d4-menu-popup [d4-name][aria-checked="true"]'))
        .map((e) => e.getAttribute('d4-name')).join(', ');
      expect(checked, 'OpenChemLib', 'the sketcher the menu checks');
      item.click();
      const switched = await host.sketcherReady();
      expect(switched.root.classList.contains('crux-sketcher'), true, 'the host\'s sketcher after Crux was picked');
      expect(grok.chem.currentSketcherType, 'Crux', 'the session\'s sketcher');
    } finally {
      dialog.close();
    }
  });

  // The SketcherBase contract in Datagrok (HOST-007 to 015, 022, 023, 042, 046, 047; spike datagrok-chem, R3): what was
  // crux-sketch's @datagrok contract feature. A user's edit is made as the user makes it, on Crux's own controls in the
  // page: a toolbar button pressed, a press on the canvas where Crux draws an atom (its `positions`).

  test('init: initialized the moment it resolves, a value set at once read at once', async () => {
    const impl = await DG.Func.find({package: 'Chem', name: 'cruxSketcher'})[0].apply() as DG.chem.SketcherBase;
    const host = new DG.chem.Sketcher();
    const box = ui.div([impl.root], {style: {width: '500px', height: '400px'}});
    const dialog = ui.dialog().add(box).show();
    try {
      await impl.init(host);
      // the moment init resolved: nothing awaited between
      const initialized = impl.isInitialized;
      impl.smiles = 'Oc1ccccc1';
      const molFile = impl.molFile;
      expect(initialized, true, 'isInitialized the moment init resolved');
      expect(canonical(molFile), canonical('Oc1ccccc1'), 'the molFile read right after the set');
      expectDrawing(impl.root.querySelector<CruxElement>('crux-sketch'), 'Oc1ccccc1');
    } finally {
      dialog.close();
      impl.detach();
    }
  });

  test('sketcherReady: on opening, and after a switch', async () => {
    grok.chem.currentSketcherType = 'Crux';
    const host = new DG.chem.Sketcher();
    host.setSmiles('CCO');
    const moment = (s: DG.chem.SketcherBase) => ({current: s === host.sketcher, initialized: s.isInitialized,
      crux: s.root.classList.contains('crux-sketcher'), type: grok.chem.currentSketcherType,
      molFile: host.getMolFile()});
    const announced: ReturnType<typeof moment>[] = [];
    host.onSketcherReady.subscribe((s) => announced.push(moment(s)));
    // asked before the host has created any implementation: it waits for the first
    const first = host.sketcherReady().then(moment);
    const dialog = ui.dialog().add(host).show();
    try {
      const opened = await first;
      expect(JSON.stringify({...opened, molFile: undefined}),
        JSON.stringify({current: true, initialized: true, crux: true, type: 'Crux'}),
        'what sketcherReady() resolved with');
      expect(canonical(opened.molFile), canonical('CCO'), 'the host\'s molecule when it resolved');
      expect(announced.length, 1, 'onSketcherReady on opening');
      // the switch the ≡ menu makes: the session's sketcher type, then the host's
      grok.chem.currentSketcherType = 'OpenChemLib';
      host.sketcherType = 'OpenChemLib';
      const switched = await host.sketcherReady().then(moment);
      expect(JSON.stringify({...switched, molFile: undefined}),
        JSON.stringify({current: true, initialized: true, crux: false, type: 'OpenChemLib'}),
        'what sketcherReady() resolved with after the switch');
      expect(canonical(switched.molFile), canonical('CCO'), 'the host\'s molecule after the switch');
      expect(announced.length, 2, 'onSketcherReady after the switch');
      expect(announced.every((a) => a.current && a.initialized && canonical(a.molFile) === canonical('CCO')), true,
        `the announcements: ${JSON.stringify(announced.map((a) => ({...a, molFile: undefined})))}`);
    } finally {
      dialog.close();
    }
  });

  for (const how of ['set', 'entered'] as const) {
    test(`a value the host takes while Crux loads is shown once ready: ${how}`, async () => {
      const value = how === 'set' ? 'CC(=O)Oc1ccccc1C(=O)O' : 'C1CCCCC1';
      grok.chem.currentSketcherType = 'Crux';
      const host = new DG.chem.Sketcher();
      let count = 0;
      let loading: boolean | null = null;
      // the host puts the new sketcher's root in its box and then awaits its init: the value comes now
      const observer = new MutationObserver(() => {
        if (host.sketcher === null || loading !== null)
          return;
        host.sketcher.onChanged.subscribe(() => count++);
        loading = !host.sketcher.isInitialized;
        if (how === 'set')
          host.setMolecule(value);
        else {
          host.molInput.value = value;
          host.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
        }
      });
      observer.observe(host.host, {childList: true});
      const dialog = ui.dialog().add(host).show();
      try {
        await host.sketcherReady();
        await delay(SETTLE_MS);
        expect(loading, true, 'Crux was loading when the value came');
        expectDrawing(host.sketcher!.root.querySelector<CruxElement>('crux-sketch'), value);
        expect(count, 1, 'onChanged events');
      } finally {
        observer.disconnect();
        dialog.close();
      }
    });
  }

  // Chem's substructure filter lost the SMILES typed into its sketcher dialog when OK closed the dialog before Crux was
  // ready (Chem's features card-state and substructure-card, run with Crux pinned; spike datagrok-chem, R2): the host
  // keeps a value set while its sketcher loads and sets it once the sketcher is ready, but a Crux detached while it
  // loaded never became ready, so nothing was set and no change event told the filter. It now answers as a sketcher
  // detached after it was ready does.
  test('a value the host takes while Crux loads is kept and said once, its dialog closed before Crux is ready',
    async () => {
    grok.chem.currentSketcherType = 'Crux';
    const host = new DG.chem.Sketcher();
    const value = 'c1ccncc1';
    let implChanges = 0;
    let hostChanges = 0;
    let moment: {loading: boolean, detached: boolean} | null = null;
    host.onChanged.subscribe(() => hostChanges++);
    const dialog = ui.dialog().add(host);
    // the host puts the new sketcher's root in its box and then awaits its init: the user types into the host's field,
    // presses Enter, and OK closes the dialog, all before Crux is ready
    const observer = new MutationObserver(() => {
      if (host.sketcher === null || moment !== null)
        return;
      const loading = host.sketcher;
      loading.onChanged.subscribe(() => implChanges++);
      host.setValue(value);
      dialog.close();
      moment = {loading: !loading.isInitialized, detached: loading.isDetached};
    });
    observer.observe(host.host, {childList: true});
    dialog.show();
    try {
      const impl = await host.sketcherReady();
      await delay(SETTLE_MS);
      expect(JSON.stringify(moment), JSON.stringify({loading: true, detached: true}),
        'Crux was loading when the value came, and detached by the dialog\'s close');
      expect(impl.root.querySelector('crux-sketch'), null, '<crux-sketch> left after the dialog closed');
      expect(implChanges, 1, 'the sketcher\'s change events');
      expect(hostChanges, 1, 'the host\'s change events');
      expect(impl.smiles, value, 'the sketcher\'s SMILES: the caller\'s string');
      expect(canonical(impl.molFile), canonical(value), 'the sketcher\'s molFile');
      expect(canonical(host.getMolFile()), canonical(value), 'the host\'s molblock');
    } finally {
      observer.disconnect();
    }
  });

  // HOST-050 failed once in Chem's Crux features (2 workers, 2026-10-07): C1CCCCC1 typed into the cell editor's field
  // and entered, then OK, and the cell kept its CCO. The Enter came while Crux loaded: the host kept the value and showed
  // it once Crux was ready, and the cell editor, on that same readiness, then gave Crux the cell's own SMILES as the
  // caller's string (`explicitMol`, so that OK without an edit leaves the cell as it was) over the one just entered.
  // Crux drew C1CCCCC1, its getters read the cell's SMILES, and OK found nothing changed.
  for (const typed of [null, 'C1CCCCC1']) {
    test(typed === null ? 'the cell editor\'s OK without an edit leaves the cell\'s SMILES as it was' :
      'the cell editor\'s OK writes a SMILES entered while Crux loads', async () => {
      // the cell editor converts the cell's SMILES with Chem's RDKit, this bundle's own module
      if (!chemCommonRdKit.moduleInitialized) {
        chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
        await chemCommonRdKit.initRdKitModuleLocal();
      }
      grok.chem.currentSketcherType = 'Crux';
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('molecule', [STEREO_SMILES, 'c1ccccc1'])]);
      const col = df.col('molecule')!;
      col.semType = DG.SEMTYPE.MOLECULE;
      col.meta.units = DG.chem.Notation.Smiles;
      // empty, so that OK adds its molecule at once (to a list it calls Chem's removeDuplicates first); put back after
      const recent = localStorage.getItem(DG.chem.Sketcher.RECENT_KEY);
      localStorage.setItem(DG.chem.Sketcher.RECENT_KEY, '[]');
      const tv = grok.shell.addTableView(df);
      let caught: (c: {host: DG.chem.Sketcher, loading: boolean}) => void;
      const created = new Promise<{host: DG.chem.Sketcher, loading: boolean}>((resolve) => caught = resolve);
      // the cell editor's Crux, as its host puts it in its box and then awaits its init: the user types into the host's
      // field and presses Enter now
      const observer = new MutationObserver((records) => {
        for (const n of records.flatMap((r) => Array.from(r.addedNodes))) {
          if (!(n instanceof HTMLElement) || !n.classList.contains('crux-sketcher'))
            continue;
          const impl = DG.Widget.find(n) as DG.chem.SketcherBase | null;
          if (!impl?.host)
            continue;
          observer.disconnect();
          const host = impl.host;
          const loading = !impl.isInitialized;
          if (typed !== null) {
            host.molInput.value = typed;
            host.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
          }
          caught({host, loading});
          return;
        }
      });
      observer.observe(document.body, {childList: true, subtree: true});
      let dialog: HTMLElement | null = null;
      try {
        await PackageFunctions.editMoleculeCell(tv.grid.cell('molecule', 0));
        const {host, loading} = await created;
        dialog = host.root.closest<HTMLElement>('.d4-dialog');
        expect(loading, true, 'Crux was loading when the cell editor opened');
        // ready: the host's and the cell editor's handlers of it have run (a microtask each), and OK comes after them
        await host.sketcherReady();
        if (typed !== null)
          expectDrawing(host.sketcher!.root.querySelector<CruxElement>('crux-sketch'), typed);
        // the OK button's handler runs after the click, and the dialog says when it is done
        const okDialog = DG.Dialog.getOpenDialogs().find((d) => d.root.contains(host.root))!;
        const okDone = new Promise<void>((resolve) => {
          const sub = okDialog.onAfterOK.subscribe(() => {
            sub.unsubscribe();
            resolve();
          });
        });
        dialog!.querySelector<HTMLElement>('[name="button-OK"]')!.click();
        await okDone;
        expect(df.get('molecule', 0), typed ?? STEREO_SMILES, 'the cell after OK');
      } finally {
        observer.disconnect();
        if (dialog !== null && document.body.contains(dialog))
          dialog.querySelector<HTMLElement>('[name="button-CANCEL"]')?.click();
        tv.close();
        if (recent === null)
          localStorage.removeItem(DG.chem.Sketcher.RECENT_KEY);
        else
          localStorage.setItem(DG.chem.Sketcher.RECENT_KEY, recent);
      }
    });
  }

  test('the getters in the change handler, in the OK handler and after the dialog closed', async () => {
    let okReads: Reads | null = null;
    const handlerReads: Reads[] = [];
    let impl: DG.chem.SketcherBase | null = null;
    const s = await openWith({onOK: () => okReads = readAll(impl!)});
    impl = s.host.sketcher!;
    impl.onChanged.subscribe(() => handlerReads.push(readAll(impl!)));
    press(s, 'bond.single');
    pressAt(s, canvasMiddle(s));
    await delay(SETTLE_MS);
    expect(handlerReads.length, 1, 'onChanged events of the bond drawn');
    expectReads(handlerReads[handlerReads.length - 1], 'CC', 'read in the change handler');
    (s.dialog.root.querySelector('[name="button-OK"]') as HTMLElement).click();
    await delay(SETTLE_MS);
    expectReads(okReads, 'CC', 'read in the OK handler');
    await awaitCheck(() => impl!.isDetached, 'Crux has not been detached', 3000);
    expect(impl.root.querySelector('crux-sketch'), null, '<crux-sketch> left after detach');
    expect(document.body.contains(impl.root), false, 'the sketcher\'s root left in the page');
    expectReads(readAll(impl), 'CC', 'read after the dialog closed');
    const smarts = await impl.getSmarts();
    expect(sameQuery(smarts!, rdkitSmarts('CC')), true, `getSmarts ${smarts} finds what CC finds`);
  });

  const setters: [string, string, string][] = [
    ['smiles', 'CC(=O)Nc1ccc(O)cc1', 'V2000'],
    ['molFile', 'CC(=O)Nc1ccc(O)cc1', 'V2000'],
    ['molV3000', 'CC(=O)Nc1ccc(O)cc1', 'V3000'],
    ['smarts', '[#6]-[#6](=[#8])-[#7,#8]', 'V2000'],
  ];
  for (const [setter, molecule, version] of setters) {
    test(`the ${setter} setter replaces the drawing, keeps the caller's string, and says so once`, async () => {
      const s = await openWith({smiles: 'c1ccccc1'});
      try {
        const {v2000, v3000} = setter === 'smarts' ? {v2000: '', v3000: ''} : molblocks(molecule);
        const text = setter === 'molFile' ? v2000 : setter === 'molV3000' ? v3000 : molecule;
        // through the host, as Datagrok sets a molecule: setMolFile hands a V3000 molblock to molV3000 by its "V3000"
        if (setter === 'smiles')
          s.host.setSmiles(text);
        else if (setter === 'smarts')
          s.host.setSmarts(text);
        else
          s.host.setMolFile(text);
        await delay(SETTLE_MS);
        expect(s.changes(), 2, 'onChanged events, the opening\'s and the setter\'s');
        const el = cruxElement(s)!;
        if (setter === 'smarts') {
          expect(el.isEmpty, false, 'Crux\'s canvas');
          expect(sameQuery(el.smarts, molecule), true, `Crux's own SMARTS ${el.smarts} finds what ${molecule} finds`);
        } else
          expectDrawing(el, molecule);
        const impl = s.host.sketcher!;
        const back = setter === 'smarts' ? await impl.getSmarts() : (impl as unknown as Record<string, string>)[setter];
        expect(back, text, `the sketcher's ${setter}: the very string the host set`);
        expect(impl.molFile.includes('V2000') && impl.molV3000.includes('V3000'), true, 'molFile V2000, molV3000 V3000');
        expect(s.host.getMolFile().includes(version), true, `the host's molblock is ${version}`);
      } finally {
        s.dialog.close();
      }
    });
  }

  const empties: [string, () => string][] = [['\'\'', () => ''], ['WHITE_MOLBLOCK', () => DG.WHITE_MOLBLOCK],
    ['WHITE_MOLBLOCK_V_3000', () => DG.WHITE_MOLBLOCK_V_3000]];
  for (const [name, empty] of empties) {
    test(`an empty value clears the canvas with one change event: ${name}`, async () => {
      const s = await openWith({smiles: 'CCO'});
      try {
        s.host.setMolecule(empty());
        await delay(SETTLE_MS);
        expect(s.changes(), 2, 'onChanged events, the opening\'s and the clear\'s');
        expectEmptyCanvas(cruxElement(s));
        expect(s.host.isEmpty(), true, 'the host reads the sketcher as empty');
      } finally {
        s.dialog.close();
      }
    });
  }

  test('an empty value clears a detached sketcher, with one change event', async () => {
    const s = await openWith({smiles: 'CCO'});
    s.dialog.close();
    await awaitCheck(() => s.host.sketcher!.isDetached, 'Crux has not been detached', 3000);
    s.host.setMolecule('');
    await delay(SETTLE_MS);
    expect(s.changes(), 2, 'onChanged events, the opening\'s and the clear\'s');
    expect(s.host.isEmpty(), true, 'the host reads the sketcher as empty');
  });

  test('the getters of a canvas the user cleared read as empty', async () => {
    const s = await openWith({smiles: 'CC(=O)O'});
    try {
      press(s, 'clear');
      await delay(SETTLE_MS);
      const impl = s.host.sketcher!;
      const empty = (m: string) => m === '' || DG.chem.Sketcher.isEmptyMolfile(m);
      expect(impl.smiles, '', 'smiles');
      expect(empty(impl.molFile), true, `molFile: ${impl.molFile}`);
      expect(empty(impl.molV3000), true, `molV3000: ${impl.molV3000}`);
      expect(await impl.getSmarts(), '', 'getSmarts');
      expect(s.host.isEmpty(), true, 'the host reads the sketcher as empty');
    } finally {
      s.dialog.close();
    }
  });

  test('the caller\'s string comes back until the user\'s first edit', async () => {
    const s = await openWith({});
    try {
      s.host.setSmiles('C1=CC=CC=C1');
      await delay(SETTLE_MS);
      expect(s.changes(), 1, 'onChanged events of the set');
      expect(s.host.sketcher!.smiles, 'C1=CC=CC=C1', 'the caller\'s string, exactly');
      press(s, 'bond.single');
      pressAt(s, atomAt(s, 0));
      await delay(SETTLE_MS);
      expect(s.changes(), 2, 'onChanged events after the user\'s bond');
      expect(s.host.sketcher!.explicitMol ?? null, null, 'the caller\'s string, after the user\'s edit');
      expect(canonical(s.host.sketcher!.smiles), canonical('Cc1ccccc1'), 'Crux\'s own SMILES');
    } finally {
      s.dialog.close();
    }
  });

  test('detach lets go of Crux without touching another, and the getters keep answering', async () => {
    const first = await openWith({smiles: 'CCO'});
    const second = await openWith({smiles: 'c1ccccc1'});
    try {
      const impl = first.host.sketcher!;
      first.dialog.close();
      await awaitCheck(() => impl.isDetached, 'the first sketcher has not been detached', 3000);
      expect(impl.root.querySelector('crux-sketch'), null, '<crux-sketch> left after detach');
      expect(document.body.contains(impl.root), false, 'the first sketcher\'s root left in the page');
      expectReads(readAll(impl), 'CCO', 'the first sketcher, read after detach');
      press(second, 'clear');
      await delay(SETTLE_MS);
      expect(second.changes(), 2, 'onChanged events of the second sketcher');
      expect(first.changes(), 1, 'onChanged events of the first sketcher');
    } finally {
      second.dialog.close();
    }
  });

  // HOST-043 (spike datagrok-platform, R2): the platform's automation finds the sketcher by its u2 name; Crux's own
  // element says it is busy until it is ready (Crux's half, its aria-busy, is crux-sketch's a11y/busy.feature)
  test('the root is named cruxSketch, and Crux\'s element is busy until it is ready', async () => {
    grok.chem.currentSketcherType = 'Crux';
    const host = new DG.chem.Sketcher();
    let busyAtInsert: string | null | undefined;
    // the moment Crux's element is put in the sketcher's root, before it is ready
    const observer = new MutationObserver(() => {
      const el = host.host.querySelector('crux-sketch');
      if (el !== null && busyAtInsert === undefined)
        busyAtInsert = el.getAttribute('aria-busy');
    });
    observer.observe(host.host, {childList: true, subtree: true});
    const dialog = ui.dialog().add(host).show();
    try {
      const sketcher = await host.sketcherReady();
      expect(sketcher.root.getAttribute('data-u2-name'), 'cruxSketch', 'the u2 name of the sketcher\'s root');
      expect(Array.from(document.querySelectorAll('[data-u2-name="cruxSketch"]')).includes(sketcher.root), true,
        'the root found by its u2 name');
      expect(busyAtInsert, 'true', 'Crux\'s aria-busy as its element was put in the page');
      expect(sketcher.root.querySelector('crux-sketch')!.getAttribute('aria-busy'), 'false',
        'Crux\'s aria-busy once the host says it is ready');
    } finally {
      observer.disconnect();
      dialog.close();
    }
  });

  // HOST-073 (R5): the platform's settles wait on the sketcher's isRenderPending and onRendered, never on a timeout
  test('onRendered fires once for each drawing Crux shows, and isRenderPending is false once it has', async () => {
    const s = await openCrux();
    try {
      const sketcher = s.host.sketcher!;
      const el = cruxElement(s)! as CruxElement & {renderCount: number};
      let rendered = 0;
      let countWhenRendered = -1;
      const sub = sketcher.onRendered.subscribe(() => {
        rendered++;
        countWhenRendered = el.renderCount;
      });
      try {
        expect(sketcher.isRenderPending, false, 'isRenderPending once Crux is ready and shown');
        const before = el.renderCount;
        s.host.setSmiles('c1ccccc1');
        await awaitCheck(() => rendered > 0, 'onRendered never fired after a write', 5000);
        await awaitCheck(() => sketcher.isRenderPending === false, 'isRenderPending stayed true', 5000);
        await delay(SETTLE_MS);
        expect(rendered, el.renderCount - before, 'onRendered against the drawings Crux showed');
        expect(rendered, 1, 'onRendered after one write');
        expect(countWhenRendered, el.renderCount, 'Crux\'s render count when onRendered fired');
        expectDrawing(el, 'c1ccccc1');
        // reading the status draws nothing
        sketcher.getWidgetStatus();
        await delay(SETTLE_MS);
        expect(rendered, 1, 'onRendered after reading the status');
      } finally {
        sub.unsubscribe();
      }
    } finally {
      s.dialog.close();
    }
  });

  test('refresh() keeps the molecule, the undo history and the caller\'s string, and says nothing', async () => {
    const s = await openWith({smiles: 'C1=CC=CC=C1'});
    const insideCanvas = (smiles: string) => {
      expectDrawing(cruxElement(s), smiles);
      const c = canvasBox(s);
      cruxElement(s)!.positions.atoms.forEach((p, i) => expect(p !== null && p.x > c.x && p.x < c.x + c.width &&
        p.y > c.y && p.y < c.y + c.height, true,
        `atom ${i} at ${JSON.stringify(p)} inside the canvas ${JSON.stringify(c)}`));
    };
    try {
      // as the substructure filter does when its panel is closed and opened again (R11)
      const root = s.host.sketcher!.root;
      root.remove();
      await new Promise((r) => requestAnimationFrame(() => r(null)));
      s.host.host.append(root);
      s.host.sketcher!.refresh();
      await delay(SETTLE_MS);
      expect(s.changes(), 1, 'onChanged events after refresh()');
      expect(s.host.sketcher!.smiles, 'C1=CC=CC=C1', 'the caller\'s string, exactly');
      insideCanvas('c1ccccc1');
      press(s, 'bond.single');
      pressAt(s, atomAt(s, 0));
      s.host.sketcher!.refresh();
      await delay(SETTLE_MS);
      expect(s.changes(), 2, 'onChanged events after the user\'s bond and refresh()');
      insideCanvas('Cc1ccccc1');
      press(s, 'undo');
      await delay(SETTLE_MS);
      expectDrawing(cruxElement(s), 'c1ccccc1');
    } finally {
      s.dialog.close();
    }
  });

  test('supportedExportFormats lists the formats Crux writes, and it writes each', async () => {
    const s = await openWith({smiles: 'CC(=O)O'});
    try {
      press(s, 'bond.single');
      pressAt(s, atomAt(s, 0));
      await delay(SETTLE_MS);
      const impl = s.host.sketcher!;
      expect(impl.supportedExportFormats.join(', '), 'smiles, mol, molV3000, smarts', 'supportedExportFormats');
      const want = canonical('CCC(=O)O');
      expect(canonical(impl.smiles), want, 'smiles');
      expect(canonical(impl.molFile) === want && impl.molFile.includes('V2000'), true, 'mol');
      expect(canonical(impl.molV3000) === want && impl.molV3000.includes('V3000'), true, 'molV3000');
      const smarts = await impl.getSmarts();
      expect(sameQuery(smarts, rdkitSmarts('CCC(=O)O')), true, `smarts ${smarts}`);
    } finally {
      s.dialog.close();
    }
  });

  test('a query\'s SMILES reads empty, with the platform\'s warning', async () => {
    const s = await openWith({smarts: SMARTS_LIST});
    try {
      expect(s.host.sketcher!.smiles, '', 'the SMILES of a query');
      await awaitCheck(() => Array.from(document.querySelectorAll('.d4-balloon.warning .d4-balloon-content'))
        .some((b) => b.textContent?.includes('Smarts cannot be converted to smiles') ?? false),
      'no warning balloon', 3000);
    } finally {
      s.dialog.close();
      DG.Balloon.closeAll();
    }
  });

  test('a malformed molblock set through the host throws nothing, and the host shows its warning', async () => {
    const validate = (molecule: string) => DG.Func.find({package: 'Chem', name: 'validateMolecule'})[0]
      .prepare({s: molecule}).callSync().getOutputParamValue();
    const s = await openWith({validation: validate, smiles: 'CCO'});
    try {
      const molblock = ['malformed', '  a bond to an atom that is not there', '',
        '  2  1  0  0  0  0  0  0  0  0999 V2000',
        '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
        '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
        '  1  5  1  0', 'M  END', ''].join('\n');
      let threw: unknown = null;
      try {
        s.host.setMolFile(molblock);
      } catch (e) {
        threw = e;
      }
      await delay(SETTLE_MS);
      expect(threw, null, 'what the host\'s call threw');
      expect(s.dialog.root.querySelector('.chem-invalid-molecule-warning')?.textContent, 'Malformed molecule',
        'the host\'s warning');
      expectDrawing(cruxElement(s), 'CCO');
      expect(s.changes(), 1, 'onChanged events');
    } finally {
      s.dialog.close();
    }
  });

  // QUERY-021 (spike datagrok-chem-hosts, R4): each query feature Crux's tools make in the substructure filter (its query
  // bonds, topology, lists, generics, the Query properties and Atom Properties' fields, a custom query), put on Crux's
  // canvas as the tool leaves it in Crux's model (read from the SMARTS, the CXSMILES label or the molfile field Crux writes
  // for it), goes the filter's own way: Crux's molblock (V2000, or V3000 with SMARTSQ groups where V2000 holds none),
  // Chem's alias translation, RDKit's query and the search. The rows it passes are the rows its intended SMARTS matches.
  // None fails, so none is turned off there (the filter's Crux turns off no query element).
  test('every query tool in the substructure filter passes the rows of its intended SMARTS', async () => {
    const probes = ['CO', 'C=O', 'CC=O', 'CCO', 'OCCO', 'c1ccncc1', 'C1CCNCC1', 'CN', 'C=N', 'CC#N', 'CCF', 'CCCl', 'CBr',
      'CI', 'C[Na]', 'C[Mg]C', 'CS', 'CC', 'C', 'C1CCCCC1', 'c1ccccc1', 'C1CC1', 'CC(C)(C)C', 'C=C', 'C#C', 'CC(=O)O',
      'Oc1ccccc1', 'C[Se]', 'CC1CCCCC1', 'Cc1ccccc1', 'NC=O', 'OC1CCCCC1', 'Oc1ccncc1'];
    const metals = '[Li,Na,K,Rb,Cs,Mg,Ca,Sr,Ba,Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg]';
    const field = (props: string, hhh = 0) => ['field', '', '', '  2  1  0  0  0  0  0  0  0  0999 V2000',
      `    0.0000    0.0000    0.0000 C   0  0  0  ${hhh}  0  0  0  0  0  0  0  0`,
      '    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0', '  1  2  1  0', ...props.split('|').filter(Boolean),
      'M  END', ''].join('\n');
    // [the tool, Crux's model of what it makes, its intended SMARTS]
    const tools: [string, string, string][] = [
      ['any bond', '[#6]~[#8]', '[#6]~[#8]'],
      ['single or double bond', '[#6]-,=[#8]', '[#6]-,=[#8]'],
      ['single or aromatic bond', '[#6]-,:[#7]', '[#6]-,:[#7]'],
      ['double or aromatic bond', '[#6]=,:[#7]', '[#6]=,:[#7]'],
      ['aromatic bond', '[#6]:[#7]', '[#6]:[#7]'],
      ['ring bond (topology)', '[#6]-&@[#6]', '[#6]-&@[#6]'],
      ['chain bond (topology)', '[#6]-&!@[#8]', '[#6]-&!@[#8]'],
      ['atom list', '[#6]-[#7,#8]', '[#6]-[#7,#8]'],
      ['NOT list', '[#6]-[!#7&!#8]', '[#6]-[!#7&!#8]'],
      ['generic A', '[#6]-[!#1]', '[#6]-[!#1]'],
      ['generic AH', 'C* |$;AH_p$|', '[#6]-*'],
      ['generic Q', 'C* |$;Q_e$|', '[#6]-[!#6&!#1]'],
      ['generic QH', 'C* |$;QH_p$|', '[#6]-[!#6]'],
      ['generic X', 'C* |$;X_p$|', '[#6]-[F,Cl,Br,I]'],
      ['generic XH', 'C* |$;XH_p$|', '[#6]-[F,Cl,Br,I,#1]'],
      ['generic M', 'C* |$;M_p$|', `[#6]-${metals}`],
      ['generic MH', 'C* |$;MH_p$|', `[#6]-${metals.replace(']', ',#1]')}`],
      ['generic *', 'C* |$;star_e$|', '[#6]-*'],
      ['ring bond count 2', field('M  RBC  1   1   2'), '[#6;x2]-[#8]'],
      ['ring bond count as drawn', field('M  RBC  1   1  -2'), '[#6;x0]-[#8]'],
      ['ring bond count 0', field('M  RBC  1   1  -1'), '[#6;x0]-[#8]'],
      ['substitution count 2', field('M  SUB  1   1   2'), '[#6;D2]-[#8]'],
      ['substitution count as drawn', field('M  SUB  1   1  -2'), '[#6;D1]-[#8]'],
      ['unsaturated', field('M  UNS  1   1   1'), '[#6;$(*=,:,#*)]-[#8]'],
      ['H count 1 or more', field('', 2), '[#6;h{1-}]-[#8]'],
      ['aromatic', '[#6;a]', '[#6;a]'],
      ['aliphatic', '[#6;A]-[#8]', '[#6;A]-[#8]'],
      ['implicit H count', '[#6;h3]-[#8]', '[#6;h3]-[#8]'],
      ['ring membership', '[#6;R1]', '[#6;R1]'],
      ['ring size', '[#6;r3]', '[#6;r3]'],
      ['connectivity', '[#6;X4]-[#6]', '[#6;X4]-[#6]'],
      ['custom query', '[#6;$([#6]=[#8])]', '[#6;$([#6]=[#8])]'],
    ];
    // the filter searches with Chem's own RDKit module, as Chem's filter tests start it
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    const col = DG.Column.fromStrings('m', probes);
    const df = DG.DataFrame.fromColumns([col]);
    col.semType = DG.SEMTYPE.MOLECULE;
    col.meta.units = DG.chem.Notation.Smiles;
    grok.chem.currentSketcherType = 'Crux';
    const filter = new SubstructureFilter();
    filter.attach(df);
    filter.applyState({columnName: 'm'});
    const dialog = ui.dialog().add(filter.root).show();
    try {
      await ui.tools.waitForElementInDom(filter.sketcher.root);
      filter.column = col;
      filter.columnName = 'm';
      filter.tableName = df.name;
      await filter.sketcher.sketcherReady();
      const impl = filter.sketcher.sketcher!;
      const el = impl.root.querySelector('crux-sketch')!;
      expect(el.mode, 'query', 'Crux\'s mode in the substructure filter');
      expect(el.disabledQueryElements.join(', '), '', 'the query elements turned off in the filter');
      const rows = (bits: DG.BitSet) => probes.filter((_, i) => bits.get(i)).join(' ');
      for (const [tool, made, intended] of tools) {
        const qmol = rdkit.get_qmol(intended);
        const want = DG.BitSet.create(probes.length, (i) => {
          const mol = rdkit.get_mol(probes[i]);
          try {
            return mol.get_substruct_match(qmol) !== '{}';
          } finally {
            mol.delete();
          }
        });
        qmol.delete();
        expect(want.trueCount > 0 && want.trueCount < probes.length, true, `${tool}: the probes tell ${intended} apart`);
        // from an empty canvas, all rows passing, so that the rows the query passes are its own
        impl.explicitMol = null;
        el.setValue('');
        await awaitCheck(() => df.filter.trueCount === probes.length, `${tool}: the filter did not clear`, 10000);
        // as the user's drawing does, the canvas's own molblock is the query (no caller's string)
        impl.explicitMol = null;
        el.setValue(made, made.includes('|') || made.includes('\n') ? undefined : {format: 'smarts'});
        expect(sameQuery(el.smarts, intended), true, `${tool}: Crux's own SMARTS ${el.smarts} finds what ${intended} finds`);
        await awaitCheck(() => rows(df.filter) === rows(want),
          `${tool}: the filter passes [${rows(df.filter)}], ${intended} matches [${rows(want)}]`, 15000);
      }
    } finally {
      dialog.close();
      filter.detach();
    }
  });
});
