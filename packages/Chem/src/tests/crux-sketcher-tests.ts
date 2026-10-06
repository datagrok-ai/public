import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

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
}

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
  grok.chem.currentSketcherType = 'Crux';
  const host = new DG.chem.Sketcher(undefined, validation);
  // As Chem's substructure filter does, before the sketcher is created
  if (substructureFilter)
    host.isSubstructureFilter = true;
  const dialog = ui.dialog().add(host).show();
  await awaitCheck(() => host.sketcher !== null, 'no sketcher was created', 15000, 10);
  let count = 0;
  host.sketcher!.onChanged.subscribe(() => count++);
  expect(host.sketcher!.root.classList.contains('crux-sketcher'), true, 'the sketcher is not Crux');
  await awaitCheck(() => host.sketcher!.isInitialized, 'Crux has not been initialized', 15000);
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
      await awaitCheck(() => host.sketcher?.isInitialized === true, 'Crux has not been initialized', 15000);
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
      await awaitCheck(() => s.host.sketcher !== crux && s.host.sketcher?.isInitialized === true,
        'OpenChemLib has not been initialized', 15000);
      await awaitCheck(() => canonical(s.host.getMolFile()) === canonical(SMILES), 'the molecule in OpenChemLib', 5000);
      const ocl = s.host.sketcher;
      grok.chem.currentSketcherType = 'Crux';
      s.host.sketcherType = 'Crux';
      await awaitCheck(() => s.host.sketcher !== ocl && s.host.sketcher?.isInitialized === true,
        'Crux has not been initialized again', 15000);
      expect(s.host.sketcher!.root.classList.contains('crux-sketcher'), true, 'the sketcher is not Crux');
      await awaitCheck(() => cruxElement(s)?.isEmpty === false, 'the molecule did not come back to Crux', 5000);
      expect(canonical(s.host.sketcher!.molFile), canonical(SMILES), 'the molecule back in Crux');
    } finally {
      s.dialog.close();
    }
  });
});
