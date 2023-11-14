import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test, expectArray, testEvent} from '@datagrok-libraries/utils/src/test';
import {
  getUserLibSettings, LibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';

category('getMolfiles', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: LibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadLibraries(true);
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true);
  });

  const csv = `Helm, atoms, bonds
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV}$$$$,11,10
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.D-Chg.dV}$$$$,11,10
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2}$$$$,9,8`;

  test('linear1', async () => {
    const df = DG.DataFrame.fromCsv(csv);
    const helmCol = df.getCol('Helm');
    const view = grok.shell.addTableView(df);
    const resCol: DG.Column<string> = await grok.functions.call('Helm:getMolfiles', {col: helmCol});

    await delay(100);
    await testEvent(view.grid.onAfterDrawContent, () => {}, () => {
      view.grid.invalidate();
    });
    expect(resCol.length, helmCol.length);
    const resMolfileList: MolfileHandlerBase[] = resCol.toList().map((molfile) => MolfileHandler.getInstance(molfile));
    expectArray(resMolfileList.map((mf) => mf.atomCount), df.getCol('atoms').toList());
    expectArray(resMolfileList.map((mf) => mf.bondCount), df.getCol('bonds').toList());
  });
});
