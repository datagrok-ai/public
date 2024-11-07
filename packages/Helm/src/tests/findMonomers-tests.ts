import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, expectObject, test} from '@datagrok-libraries/utils/src/test';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings,
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {findMonomers, parseHelm} from '../utils';
import {initHelmMainPackage} from './utils';


/** Tests with default monomer library */
category('findMonomers', () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    await initHelmMainPackage();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Tests 'findMonomers' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries for tests
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
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
    const monomerLib = monomerLibHelper.getMonomerLib();
    const monomerSymbolList: string[] = parseHelm(testHelmValue);
    const resMissedSet: Set<string> = findMonomers(monomerSymbolList, monomerLib);
    expectObject(resMissedSet, tgtMissedSet);
  }
});
