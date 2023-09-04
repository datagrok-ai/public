import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, delay, expect} from '@datagrok-libraries/utils/src/test';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import { KetcherSketcher } from '../ketcher';

const testSmiles = 'c1ccccc1';
const testMolfile = `
  Ketcher  8 92217 42D 1   1.00000     0.00000     0                      

  6  5  0  0  0  0  0  0  0  0999 V2000
   -3.6828    1.5285    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3973    1.1160    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1117    1.5285    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6828    2.3535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9683    1.1160    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2339    1.4893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  2  0  0  0  0
  1  5  1  0  0  0  0
  5  6  1  0  0  0  0
M  END`;
const testSmarts = '[!#6&!#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';

export async function _testSetSmiles() {
    const sketcher = await createKetcher();
    sketcher.setKetcherMolecule(testSmiles);
    await awaitCheck(() => sketcher._smiles === testSmiles, 'molecule has not been set', 5000);
}

export async function _testSetMolfile() {
    const sketcher = await createKetcher();
    sketcher.setKetcherMolecule(testMolfile);
    await awaitCheck(() => !!sketcher._molV2000 && (sketcher._molV2000!.split('\n').splice(3).join('') === testMolfile.split('\n').splice(3).join('')), 
        'molecule has not been set', 5000);
}

export async function _testSetSmarts() {
    const sketcher = await createKetcher();
    const t = new Promise((resolve, reject) => {
        sketcher!.onChanged.subscribe(async (_: any) => {
            try {
                const resultMol = await sketcher.getSmarts();
                resolve(resultMol);
            } catch (error) {
                reject(error);
            }
        });
    });
    sketcher.setKetcherMolecule(testMolfile);
    const res = await t;
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const qmolTest = rdkit.get_qmol(testSmarts);
    const smarts = qmolTest.get_smarts();
    qmolTest?.delete();
    const qmol = rdkit.get_qmol(res);
    const resSmarts = qmol.get_smarts();
    qmol?.delete();
    expect(resSmarts, smarts);
}

async function createKetcher(): Promise<KetcherSketcher> {
    const sketcher = await DG.Func.find({tags: ['moleculeSketcher'], name: 'ketcherSketcher'})[0].apply();
    const sketcherHost = new Sketcher();
    await sketcher.init(sketcherHost);
    await delay(1000);
    return sketcher;
}