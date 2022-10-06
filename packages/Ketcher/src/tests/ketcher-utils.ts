import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;

const testSmiles = 'c1ccccc1';
const testMolfile = `
  Ketcher  8 92217 42D 1   1.00000     0.00000     0                      

  6  5  0  0  0  0            999 V2000
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
M  END
`;
const testSmarts = '[!#6&!#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';

export async function _testSetSmiles() {
    const sketcher = await createKetcherAndSetMolecule(testSmiles);
    const smilesRes = await sketcher._sketcher!.getSmiles();
    expect(smilesRes, testSmiles);
}

export async function _testSetMolfile() {
    const sketcher = await createKetcherAndSetMolecule(testMolfile);
    const res = await sketcher._sketcher!.getMolfile();
    expect(res.split('\n').splice(3).join(''), testMolfile.split('\n').splice(3).join(''));
}

export async function _testSetSmarts() {
    const sketcher = await createKetcherAndSetMolecule(testSmarts);
    const res = await sketcher._sketcher!.getSmarts();
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const qmolTest = rdkit.get_qmol(testSmarts);
    const smarts = qmolTest.get_smarts();
    qmolTest?.delete();
    const qmol = rdkit.get_qmol(res);
    const resSmarts = qmol.get_smarts();
    qmol?.delete();
    expect(resSmarts, smarts);
}

async function createKetcherAndSetMolecule(molecule: string) {
    const sketcher = await DG.Func.find({tags: ['moleculeSketcher'], name: 'ketcherSketcher'})[0].apply();
    const sketcherHost = new Sketcher();
    await sketcher.init(sketcherHost);
    await delay(1000);
    await sketcher._sketcher!.setMolecule(molecule);
    return sketcher;
}