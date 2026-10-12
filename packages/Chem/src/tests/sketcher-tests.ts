import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import {category, expect, test, before, after, testEvent, delay, awaitCheck} from '@datagrok-libraries/test/src/test';
import {malformedMolblock, molV2000, molV3000} from './utils';


category('sketcher testing', () => {
  let rdkitModule: any;
  let funcs: DG.Func[];

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
    funcs = DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}});
    await sketchersWarmUp(funcs);
    grok.shell.closeAll();
  });

  test('smiles', async () => {
    await testSmiles(rdkitModule, funcs);
  });

  test('input_smiles', async () => {
    await testSmiles(rdkitModule, funcs, true);
  });

  test('molV2000', async () => {
    await testMolblock(rdkitModule, funcs, 'V2000');
  });

  test('paste_input_molV2000', async () => {
    await testMolblock(rdkitModule, funcs, 'V2000', true);
  });

  test('molV3000', async () => {
    await testMolblock(rdkitModule, funcs, 'V3000');
  });

  test('paste_input_molV3000', async () => {
    await testMolblock(rdkitModule, funcs, 'V3000', true);
  });

  test('smarts', async () => {
    await testSmarts(rdkitModule, funcs);
  });

  test('inchi', async () => {
    await testInchi(rdkitModule, funcs);
  }, {timeout: 90000});

  test('malformed input', async () => {
    await testMolblock(rdkitModule, funcs, 'V2000', false, true);
  });

  // With an empty browser cache, applying Ketcher's function loads its whole package, for seconds; a sketcher chosen
  // meanwhile (the options menu) is shown first, and Ketcher, when it arrived, replaced it in the host.
  test('a sketcher chosen while Ketcher loads stays', async () => {
    const ketcher = funcs.find((f) => f.name === 'ketcherSketcher');
    if (!ketcher)
      return;
    let release = () => {};
    const held = new Promise<void>((resolve) => release = resolve);
    let applying = false;
    let applied: Promise<chem.SketcherBase> | null = null;
    const [s, d] = await sketcherSwitchingFrom(ketcher, (f) => async () => {
      applying = true;
      await held;
      applied = f.apply();
      return applied;
    });
    await awaitCheck(() => applying, 'the function of Ketcher has not been applied', 10000);
    const chosen = await switchSketcher(s, 'OpenChemLib');
    release();
    await applied;
    await delay(500);
    expect(s.sketcher === chosen, true);
    expect(s.host.contains(chosen.root), true);
    expect(s.host.querySelector('.ketcher-host') === null, true);
    d.close();
  });

  // The same while Ketcher initializes (its editor's first load): once it was ready, it took the host's change
  // subscription from the sketcher chosen meanwhile, and set the molecule it started with over that one's.
  test('a sketcher chosen while Ketcher initializes keeps its molecule', async () => {
    const ketcher = funcs.find((f) => f.name === 'ketcherSketcher');
    if (!ketcher)
      return;
    let release = () => {};
    const held = new Promise<void>((resolve) => release = resolve);
    let initializing = false;
    let initialized: Promise<void> | null = null;
    const [s, d] = await sketcherSwitchingFrom(ketcher, (f) => async () => {
      const implementation: chem.SketcherBase = await f.apply();
      const init = implementation.init.bind(implementation);
      implementation.init = async (host: chem.Sketcher) => {
        initializing = true;
        await held;
        initialized = init(host);
        return initialized;
      };
      return implementation;
    });
    await awaitCheck(() => initializing, 'Ketcher has not begun to initialize', 10000);
    const chosen = await switchSketcher(s, 'OpenChemLib');
    let changes = 0;
    s.onChanged.subscribe(() => changes++);
    s.setSmiles(exampleInchiSmiles);
    await awaitCheck(() => compareTwoMols(rdkitModule, rdkitModule.get_mol(exampleInchiSmiles), s.getMolFile()),
      'the chosen sketcher does not show its molecule', 5000);
    release();
    await initialized;
    await delay(500);
    expect(s.sketcher === chosen, true);
    expect(s.host.contains(chosen.root), true);
    expect(compareTwoMols(rdkitModule, rdkitModule.get_mol(exampleInchiSmiles), s.getMolFile()), true);
    const before = changes;
    s.setSmiles(exampleSmiles);
    await awaitCheck(() => changes > before, 'the host hears no change of the chosen sketcher', 5000);
    await delay(300);
    expect(changes - before, 1);
    d.close();
  });
});

/** A host in a dialog, switching to `from` (as the session's sketcher), whose function's apply is `apply(from)`. */
async function sketcherSwitchingFrom(from: DG.Func, apply: (f: DG.Func) => () => Promise<chem.SketcherBase>):
  Promise<[Sketcher, DG.Dialog]> {
  chem.currentSketcherType = from.friendlyName;
  const s = new Sketcher();
  s.createSketcher();
  // the switch looks its function up after a turn: replaced here, before it does
  s.sketcherFunctions = s.sketcherFunctions.map((f) => f.name !== from.name ? f :
    ({name: f.name, friendlyName: f.friendlyName, apply: apply(f)} as unknown as DG.Func));
  const d = ui.dialog().add(s).show();
  return [s, d];
}

/** Chooses `friendlyName` as the options menu does, and resolves to its implementation once ready. */
async function switchSketcher(s: Sketcher, friendlyName: string): Promise<chem.SketcherBase> {
  chem.currentSketcherType = friendlyName;
  s.sketcherType = friendlyName;
  return await s.sketcherReady();
}


function compareTwoMols(rdkitModule: any, mol: any, resMolfile: any): boolean {
  const mol2 = rdkitModule.get_mol(resMolfile);
  const match1 = mol.get_substruct_match(mol2);
  const match2 = mol2.get_substruct_match(mol);
  mol2.delete();
  return match1 !== '{}' && match2 !== '{}';
}

async function testSmarts(rdkitModule: any, funcs: DG.Func[]) {
  const mol = rdkitModule.get_mol(exampleSmiles);
  const qmol = rdkitModule.get_qmol(convertedSmarts);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher') continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await s.sketcherReady();
    const t = new Promise((resolve, reject) => {
      s.sketcher!.onChanged.subscribe(async (_: any) => {
        try {
          const resultSmarts = await s.getSmarts();
          resolve(resultSmarts);
        } catch (error) {
          reject(error);
        }
      });
    });
    setTimeout(()=> {s.setSmarts(convertedSmarts);}, 1000);
    const resSmarts = await t;
    const qmol2 = rdkitModule.get_qmol(resSmarts);
    const match1 = mol.get_substruct_match(qmol2);
    const match2 = mol.get_substruct_match(qmol);
    expect(match1 === match2, true);
    qmol2.delete();
    d.close();
  }
  mol?.delete();
  qmol?.delete();
}

const validationFunc = (s: string) => {
  const valFunc = DG.Func.find({package: 'Chem', name: 'validateMolecule'})[0];
  const funcCall: DG.FuncCall = valFunc.prepare({s});
  funcCall.callSync();
  const res = funcCall.getOutputParamValue();
  return res;
};

async function testSmiles(rdkitModule: any, funcs: DG.Func[], input?: boolean, malformed?: boolean) {
  const mol = rdkitModule.get_mol(exampleSmiles);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher(undefined, validationFunc);
    const d = ui.dialog().add(s).show();
    await s.sketcherReady();
    if (input) {
      setTimeout(() => {
        s.molInput.value = exampleSmiles;
        s.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
      }, 1000);
    } else
      setTimeout(() => s.setSmiles(exampleSmiles), 1000);
    if (malformed) {
      await awaitCheck(() => {
        const elements = d.root.getElementsByClassName('chem-invalid-molecule-warning');
        return elements.length > 0 && elements[0].children.length > 0 && elements[0].children[0].textContent === 'Malformed molecule';
      }, 'error div has not been created', 10000, 500);
    } else {
      await awaitCheck(() => {
        const resMolblock = input ? s.getMolFile() : s.getSmiles();
        return compareTwoMols(rdkitModule, mol, resMolblock);
      }, 'mols are not equal', 3000);
    }
    d.close();
  }
  mol?.delete();
}

async function testMolblock(rdkitModule: any, funcs: DG.Func[], ver: string, input?: boolean, malformed?: boolean) {
  const molfile = ver === 'V2000' ? molV2000 : molV3000;
  const mol = rdkitModule.get_mol(molfile);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher(undefined, validationFunc);
    const d = ui.dialog().add(s).show();
    await s.sketcherReady();
    if (input) {
      setTimeout(() => {
        let dT = null;
        try {dT = new DataTransfer();} catch (e) { }
        const evt = new ClipboardEvent('paste', {clipboardData: dT});
        evt.clipboardData!.setData('text/plain', molfile);
        s.molInput.value = molfile;
        s.molInput.dispatchEvent(evt);
      }, 1000);
    } else
      setTimeout(() => s.setMolFile(malformed ? malformedMolblock : molfile), 1000);
    if (malformed) {
      await awaitCheck(() => {
        const elements = d.root.getElementsByClassName('chem-invalid-molecule-warning');
        return elements.length > 0 && elements[0].children.length > 0 && elements[0].children[0].textContent === 'Malformed molecule';
      }, 'error div has not been created', 10000, 500);
    } else {
      await awaitCheck(() => {
        const resMolblock = s.getMolFile();
        return compareTwoMols(rdkitModule, mol, resMolblock);
      }, 'mols are not equal', 3000);
    }
    d.close();
  }
  mol?.delete();
}

async function testInchi(rdkitModule: any, funcs: DG.Func[]) {
  const mol = rdkitModule.get_mol(exampleInchiSmiles);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await s.sketcherReady();
    s.molInput.value = exampleInchi;
    s.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
    await awaitCheck(() => {
      const resMolblock = s.getMolFile();
      return compareTwoMols(rdkitModule, mol, resMolblock);
    }, 'mols are not equal', 3000);
    d.close();
  }
  mol?.delete();
}

export async function sketchersWarmUp(funcs: DG.Func[]) {
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await s.sketcherReady();
  }
  chem.currentSketcherType = 'OpenChemLib';
}

const exampleSmiles = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const convertedSmarts = '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
const exampleInchi = 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H';
const exampleInchiSmiles = 'c1ccccc1';
