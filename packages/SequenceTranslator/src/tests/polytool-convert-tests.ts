import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {doPolyToolConvert} from '../polytool/pt-conversion';
import {getRules} from '../polytool/pt-rules';

import {_package} from '../package-test';


category('PolyTool: Convert', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; //backup

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
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
        'R-F-D(2)-T-G-H-F-Y-P-NH2(2)',
      ],
      tgt: [
        'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$',
        'PEPTIDE1{C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$',
        'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$',
        'PEPTIDE1{C.T.G.H.F.H.P.C}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$',
        'PEPTIDE1{R.F.D.T.G.H.F.Y.P.[NH2]}$PEPTIDE1,PEPTIDE1,10:R2-3:R3$$$',
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


  test('ui-col-wo-table', async () => {
    const packageName = _package.name;
    const seqCol = DG.Column.fromStrings('seq', tests['cyclized'].src);
    const helmCol: DG.Column = await grok.functions.call(`${packageName}:polyToolConvert2`, {
      seqCol: seqCol,
      generateHelm: true,
      chiralityEngine: true,
      rules: ['rules_example.json']
    });
    expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helmCol.meta.units, NOTATION.HELM);
    expect(helmCol.dataFrame, null);
  });

  test('ui-col', async () => {
    const packageName = _package.name;
    const seqCol = DG.Column.fromStrings('seq', tests['cyclized'].src);
    const df = DG.DataFrame.fromColumns([seqCol]);
    const tv = grok.shell.addTableView(df);

    const helmCol: DG.Column = await grok.functions.call(`${packageName}:polyToolConvert2`, {
      seqCol: seqCol,
      generateHelm: true,
      chiralityEngine: true,
      rules: ['rules_example.json']
    });
    expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helmCol.meta.units, NOTATION.HELM);
    expect(helmCol.dataFrame.id, df.id);

    await testEvent(tv.grid.onAfterDrawContent, () => {}, async () => {
      tv.grid.invalidate();
    }, 15000);
    expect(helmCol.getTag(DG.TAGS.CELL_RENDERER), 'helm');
  });
});
