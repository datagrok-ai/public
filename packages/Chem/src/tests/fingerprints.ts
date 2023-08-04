import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, category, test} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { chemGetFingerprints } from '../chem-searches';
import { Fingerprint } from '../utils/chem-common';

category('fingerprints', async () => {
  let molecules: DG.DataFrame;
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molecules = DG.Test.isInBenchmark ? await grok.data.files.openTable("Demo:Files/chem/smiles_1M.zip") :
      grok.data.demo.molecules(100);
  });

  test('Pattern', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Pattern, false)!;
  });

  test('Morgan', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Morgan, false)!;
  });

});
