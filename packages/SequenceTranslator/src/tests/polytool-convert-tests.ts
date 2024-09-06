import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {doPolyToolConvert} from '../polytool/pt-conversion';
import {getRules} from '../polytool/pt-rules';


category('PolyTool: Convert', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; //backup

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    // Clear settings to test default
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadMonomerLib(true);
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const tests: { [testName: string]: { src: string[], tgt: string[] } } = {
    'cyclized': {
      src: [
        'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI',
        'C(1)-T-G-H-F-Y-P-C(1)-meI',
        'R-F-C(1)-T-G-H-F-Y-P-C(1)',
        'C(1)-T-G-H-F-H-P-C(1)',
        'R-F-L(2)-T-G-H-F-Y-P-NH2(2)',
      ],
      tgt: [
        'PEPTIDE1{[R].[F].[C].[T].[G].[H].[F].[Y].[P].[C].[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$',
        'PEPTIDE1{[C].[T].[G].[H].[F].[Y].[P].[C].[meI]}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$',
        'PEPTIDE1{[R].[F].[C].[T].[G].[H].[F].[Y].[P].[C]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$',
        'PEPTIDE1{[C].[T].[G].[H].[F].[H].[P].[C]}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$',
        'PEPTIDE1{[R].[F].[L].[T].[G].[H].[F].[Y].[P].[NH2]}$PEPTIDE1,PEPTIDE1,10:R2-3:R3$$$',
      ]
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const res = doPolyToolConvert(testData.src, rules);
      expectArray(res, testData.tgt);
    });
  }
});
