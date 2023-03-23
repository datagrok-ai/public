import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import {category, expect, test, before, after, testEvent, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {molV2000, molV3000} from './utils';


category('sketcher testing', () => {
  let rdkitModule: any;
  let funcs: DG.Func[];

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
    funcs = DG.Func.find({tags: ['moleculeSketcher']});
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
  }, {skipReason: 'GROK-12588'});

  after(async () => {
    grok.shell.closeAll();
  });
});


function compareTwoMols(rdkitModule: any, mol: any, resMolfile: any) {
  const mol2 = rdkitModule.get_mol(resMolfile);
  const match1 = mol.get_substruct_match(mol2);
  expect(match1 !== '{}', true);
  const match2 = mol2.get_substruct_match(mol);
  expect(match2 !== '{}', true);
  mol2.delete();
}

async function testSmarts(rdkitModule: any, funcs: DG.Func[]) {
  const mol = rdkitModule.get_mol(exampleSmiles);
  const qmol = rdkitModule.get_qmol(convertedSmarts);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher') continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await awaitCheck(() => s.sketcher !== null, undefined, 5000);
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

async function testSmiles(rdkitModule: any, funcs: DG.Func[], input?: boolean) {
  const mol = rdkitModule.get_mol(exampleSmiles);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await awaitCheck(() => s.sketcher !== null, undefined, 5000);
    const t = new Promise((resolve, reject) => {
      s.sketcher!.onChanged.subscribe(async (_: any) => {
        try {
          const resultMol = input ? s.getMolFile() : s.getSmiles();
          resolve(resultMol);
        } catch (error) {
          reject(error);
        }
      });
    });
    if (input) {
      setTimeout(() => {
        s.molInput.value = exampleSmiles;
        s.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
      }, 1000);
    } else setTimeout(() => s.setSmiles(exampleSmiles), 1000);
    const resMolblock = await t;
    compareTwoMols(rdkitModule, mol, resMolblock);
    d.close();
  }
  mol?.delete();
}

async function testMolblock(rdkitModule: any, funcs: DG.Func[], ver: string, input?: boolean) {
  const molfile = ver === 'V2000' ? molV2000 : molV3000;
  const mol = rdkitModule.get_mol(molfile);
  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await awaitCheck(() => s.sketcher !== null, undefined, 5000);
    const t = new Promise((resolve, reject) => {
      s.sketcher!.onChanged.subscribe(async (_: any) => {
        try {
          const resultMol = s.getMolFile();
          resolve(resultMol);
        } catch (error) {
          reject(error);
        }
      });
    });
    if (input) {
      setTimeout(() => {
        let dT = null;
        try {dT = new DataTransfer();} catch (e) { }
        const evt = new ClipboardEvent('paste', {clipboardData: dT});
        evt.clipboardData!.setData('text/plain', molfile);
        s.molInput.value = molfile;
        s.molInput.dispatchEvent(evt);
      }, 1000);
    } else setTimeout(() => s.setMolFile(molfile), 1000);
    const resMolblock = await t;
    compareTwoMols(rdkitModule, mol, resMolblock);
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
    await awaitCheck(() => s.sketcher !== null, undefined, 5000);
    s.molInput.value = exampleInchi;
    await delay(5000);
    await testEvent(s.onChanged,
      () => {compareTwoMols(rdkitModule, mol, s.getMolFile());},
      () => {s.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));}, 10000);
    d.close();
  }
  mol?.delete();
}

const exampleSmiles = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const convertedSmarts = '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
const exampleInchi = 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H';
const exampleInchiSmiles = 'c1ccccc1';
