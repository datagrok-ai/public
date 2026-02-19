import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {awaitCheck, expect} from '@datagrok-libraries/test/src/test';
import Sketcher = grok.chem.Sketcher;
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
    const {sketcher, dialog} = await createKetcher();
    sketcher.setKetcherMolecule(testSmiles);
    await awaitCheck(() => sketcher._smiles === testSmiles, 'molecule has not been set', 10000);
    dialog.close();
}

export async function _testSetMolfile() {
    const {sketcher, dialog} = await createKetcher();
    sketcher.setKetcherMolecule(testMolfile);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    await awaitCheck(() => {
        if (!!sketcher._molV2000) {
            const mol1 = rdkit.get_mol(testMolfile);
            const mol2 = rdkit.get_mol(sketcher._molV2000);
            const match1 = mol1.get_substruct_match(mol2);
            const match2 = mol2.get_substruct_match(mol1);
            mol1.delete();
            mol2.delete();
            return match1 !== '{}' && match2 !== '{}';
        }
        return false;
    }, 
        'molecule has not been set', 10000);
    dialog.close();
}

export async function _testSetSmarts() {
    const {sketcher, dialog} = await createKetcher();
    sketcher.setKetcherMolecule(testSmarts);
    await awaitCheck(() => sketcher._molV3000 !== null && sketcher._molV3000 !== '', 'molecule has not been set', 10000);
    const res = await sketcher.getSmarts();
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const qmolTest = rdkit.get_qmol(testSmarts);
    const smarts = qmolTest.get_smarts();
    qmolTest?.delete();
    const qmol = rdkit.get_qmol(res);
    const resSmarts = qmol.get_smarts();
    qmol?.delete();
    expect(resSmarts, smarts);
    dialog.close();
}

async function createKetcher(): Promise<any> {
    const func = await DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}, name: 'ketcherSketcher'})[0];
    grok.chem.currentSketcherType = func.friendlyName;
    const s = new Sketcher();
    const d = ui.dialog().add(s).show();
    await awaitCheck(() => s.sketcher?.isInitialized === true, undefined, 10000);
    return {sketcher: s.sketcher, dialog: d};
}
