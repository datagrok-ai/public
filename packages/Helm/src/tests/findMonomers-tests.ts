import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package-test';
import {after, before, category, delay, expect, expectObject, test} from '@datagrok-libraries/utils/src/test';
import {findMonomers} from '../utils';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

const LIB_STORAGE_NAME = 'Libraries';
export const LIB_DEFAULT: { [fileName: string]: string } = {'HELMCoreLibrary.json': 'HELMCoreLibrary.json'};


/** Tests with default monomer library */
category('findMonomers', () => {

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibrariesSettings: any = null;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibrariesSettings = await grok.dapi.userDataStorage.get(LIB_STORAGE_NAME, true);

    // Tests 'findMonomers' requires default monomer library loaded
    await grok.dapi.userDataStorage.post(LIB_STORAGE_NAME, LIB_DEFAULT, true);
    await monomerLibHelper.loadLibraries(true); // load default libraries
  });

  after(async () => {
    await grok.dapi.userDataStorage.put(LIB_STORAGE_NAME, userLibrariesSettings, true);
  });

  const tests: { [testName: string]: { test: string, tgt: Set<string> } } = {
    'withoutMissed': {
      test: 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$',
      tgt: new Set<string>(),
    },
    'withMissed':
      {
        test: 'PEPTIDE1{meI.missed2.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$',
        tgt: new Set<string>(['missed2'])
      }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(testName, async () => {
      _testFindMonomers(testData.test, testData.tgt);
    });
  }

  function _testFindMonomers(testHelmValue: string, tgtMissedSet: Set<string>): void {
    const resMissedSet: Set<string> = findMonomers(testHelmValue);
    expectObject(resMissedSet, tgtMissedSet);
  }
});
