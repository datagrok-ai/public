import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test, expectArray, testEvent} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {
  getUserLibSettings, setUserLibSettings,
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileHandlerBase} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler-base';

import {awaitGrid, initHelmMainPackage} from './utils';

import {_package} from '../package-test';

category('getMolfiles', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    await initHelmMainPackage();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const csv = `Helm, atoms, bonds
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV}$$$$,11,10
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.D-Chg.dV}$$$$,11,10
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2}$$$$,9,8`;

  test('linear1', async () => {
    const df = DG.DataFrame.fromCsv(csv);
    const helmCol = df.getCol('Helm');
    const view = grok.shell.addTableView(df);

    _package.logger.debug('Helm/getMolfiles/linear1, grid awaiting 1');
    await awaitGrid(view.grid, 5000);
    _package.logger.debug('Helm/getMolfiles/linear1, grid awaited 1');

    const resCol: DG.Column<string> = await grok.functions.call('Helm:getMolfiles', {col: helmCol});
    expect(resCol.length, helmCol.length);
    const resMolfileList: MolfileHandlerBase[] = resCol.toList().map((molfile) => MolfileHandler.getInstance(molfile));
    expectArray(resMolfileList.map((mf) => mf.atomCount), df.getCol('atoms').toList());
    expectArray(resMolfileList.map((mf) => mf.bondCount), df.getCol('bonds').toList());

    _package.logger.debug('Helm/getMolfiles/linear1, grid awaiting 2');
    await awaitGrid(view.grid, 5000);
    _package.logger.debug('Helm/getMolfiles/linear1, grid awaited 2');
  });
});
