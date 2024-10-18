import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {chemGetFingerprints} from '../chem-searches';
import {Fingerprint} from '../utils/chem-common';

category('fingerprints', async () => {
  let molecules: DG.DataFrame;
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molecules = DG.Test.isInBenchmark ? await grok.data.files.openTable('System:AppData/Chem/tests/smiles_100K.zip') :
      grok.data.demo.molecules(100);
  });

  test('Pattern', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Pattern, false)!;
  }, {benchmark: true});

  test('Morgan', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Morgan, false)!;
  }, {benchmark: true});

  test('fp_withoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    const res = await chemGetFingerprints(col, Fingerprint.Morgan, false)!;
    expect(res.length, 2);
  }, {benchmark: true});
});
