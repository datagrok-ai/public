import {category, expect, test, before, after} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import {_package} from '../package-test';
const exampleSmiles = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const convertedSmarts = '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
const exampleInchi = 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H';
const exampleInchiSmiles = 'c1ccccc1';

category('sketcher testing', () => {

  let rdkitModule: any;
  let funcs: DG.Func[];

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
    funcs = DG.Func.find({tags: ['moleculeSketcher']});
  });
  
  test('smiles', async () => {
    await testSmiles(rdkitModule, funcs);
  });

  test('input_smiles', async () => {
    await testSmiles(rdkitModule, funcs, true);
  });

  test('molV2000', async () => {
    await testMolV2000(rdkitModule, funcs);
  });

  test('paste_input_molV2000', async () => {
    await testMolV2000(rdkitModule, funcs, true);
  });

  test('smarts', async () => {
    await testSmarts(rdkitModule, funcs);
  });

  test('inchi', async () => {
    await testInchi(rdkitModule, funcs);
  });

  after(async () => {
    
  });

});

async function initSketcher(sw: Sketcher) {
  const t = new Promise(async (resolve, reject) => {
    sw.sketcherCreated.subscribe(async (_: any) => {
      try {
        resolve(true);
      } catch (error) {
        reject(error);
      }
    });
  });
  await t;
}

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
    if (func.name === 'chemDrawSketcher')
      continue;
    await window.localStorage.setItem('sketcher', func.name);
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    if (!s.sketcher)
      await initSketcher(s);
    setTimeout(()=> {s.setSmarts(convertedSmarts)}, 1000);
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
    await window.localStorage.setItem('sketcher', func.name);
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    if (!s.sketcher)
      await initSketcher(s);
    if (input) {
      s.molInput.value = exampleSmiles;
      s.molInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
    } else {
      setTimeout(() => {s.setSmiles(exampleSmiles)}, 1000);
    }
    const t = new Promise((resolve, reject) => {
      s.sketcher!.onChanged.subscribe(async (_: any) => {
        try {
          const resultMol = input ? s.getMolFile(): s.getSmiles();
          resolve(resultMol);
        } catch (error) {
          reject(error);
        }
      });
    });
    const resMolblock = await t;
    compareTwoMols(rdkitModule, mol, resMolblock);
    d.close();
  }
  mol?.delete();
}

async function testMolV2000(rdkitModule: any, funcs: DG.Func[], input?: boolean) {
  const data = DG.DataFrame.fromCsv(await _package.files.readAsText('test.csv'));
  const molfileV2000 = data.get('molecule', 0);
  const mol = rdkitModule.get_mol(molfileV2000);

  for (const func of funcs) {
    if (func.name === 'chemDrawSketcher')
      continue;
    await window.localStorage.setItem('sketcher', func.name);
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    if (!s.sketcher)
      await initSketcher(s);
    if (input) {
      setTimeout(() => {
        let dT = null;
        try { dT = new DataTransfer(); } catch (e) { }
        var evt = new ClipboardEvent('paste', { clipboardData: dT });
        evt.clipboardData!.setData('text/plain', molfileV2000);
        s.molInput.value = molfileV2000;
        s.molInput.dispatchEvent(evt);
      }, 1000);
    } else {
      setTimeout(() => {s.setMolFile(molfileV2000)}, 1000);
    }
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
    await window.localStorage.setItem('sketcher', func.name);
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    if (!s.sketcher)
      await initSketcher(s);
    s.molInput.value = exampleInchi;
    s.molInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
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
    const resMolblock = await t;
    compareTwoMols(rdkitModule, mol, resMolblock);
    d.close();
  }
  mol?.delete();
}
