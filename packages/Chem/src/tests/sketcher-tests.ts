import {category, expect, test, before, after, testEvent, delay} from '@datagrok-libraries/utils/src/test';
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
const exampleMol = `
Accelrys05311914342D 1   1.00000     0.00000     0

 18 19  0     0  0            999 V2000
    2.9291   -5.8667    0.0000 C   0  0  2  0  0  0           0  0  0
    3.7541   -5.8667    0.0000 C   0  0  2  0  0  0           0  0  0
    4.0109   -5.0826    0.0000 O   0  0  0  0  0  0           0  0  0
    3.3416   -4.5958    0.0000 C   0  0  2  0  0  0           0  0  0
    2.6766   -5.0826    0.0000 C   0  0  0  0  0  0           0  0  0
    3.3404   -3.7708    0.0000 N   0  0  3  0  0  0           0  0  0
    4.2383   -6.5347    0.0000 O   0  0  0  0  0  0           0  0  0
    2.4433   -6.5335    0.0000 N   0  0  0  0  0  0           0  0  0
    1.6229   -6.4464    0.0000 N   0  3  0  0  0  0           0  0  0
    5.0589   -6.4494    0.0000 C   0  0  0  0  0  0           0  0  0
    0.7983   -6.3826    0.0000 N   0  5  0  0  0  0           0  0  0
    4.0576   -3.3612    0.0000 C   0  0  0  0  0  0           0  0  0
    4.0583   -2.5398    0.0000 C   0  0  0  0  0  0           0  0  0
    3.3451   -2.1245    0.0000 C   0  0  0  0  0  0           0  0  0
    2.6294   -2.5369    0.0000 N   0  0  0  0  0  0           0  0  0
    2.6270   -3.3645    0.0000 C   0  0  0  0  0  0           0  0  0
    3.3469   -1.2995    0.0000 O   0  0  0  0  0  0           0  0  0
    1.9131   -3.7781    0.0000 O   0  0  0  0  0  0           0  0  0
  8  9  2  0     0  0
  4  5  1  0     0  0
  7 10  1  0     0  0
  5  1  1  0     0  0
  9 11  2  0     0  0
  6 12  1  0     0  0
  1  2  1  0     0  0
  4  6  1  6     0  0
  2  7  1  6     0  0
  2  3  1  0     0  0
  6 16  1  0     0  0
 12 13  2  0     0  0
 13 14  1  0     0  0
 14 15  1  0     0  0
 15 16  1  0     0  0
  1  8  1  1     0  0
 14 17  2  0     0  0
  3  4  1  0     0  0
 16 18  2  0     0  0
M  CHG  2   9   1  11  -1
M  END`

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
  const molfileV2000 = exampleMol;
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
    //s.molInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
    // const t = new Promise((resolve, reject) => {
    //   s.sketcher!.onChanged.subscribe(async (_: any) => {
    //     try {
    //       const resultMol = s.getMolFile();
    //       resolve(resultMol);
    //     } catch (error) {
    //       reject(error);
    //     }
    //   });
    // });
    await delay(5000);
    await testEvent(s.sketcher!.onChanged,
                    () => {compareTwoMols(rdkitModule, mol, s.getMolFile())}, 
                    () => {s.molInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }))}, 5000);
    // const resMolblock = await t;
    // compareTwoMols(rdkitModule, mol, resMolblock);
    d.close();
  }
  mol?.delete();
}
