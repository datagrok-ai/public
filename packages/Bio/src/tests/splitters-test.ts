import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';
import {TAGS as bioTAGS, splitterAsFasta} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {splitterAsHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {splitToMonomersUI} from '../utils/split-to-monomers';
import {awaitGrid} from './utils';
import * as C from '../utils/constants';
import {getHelmMonomers} from '../package';

import {_package} from '../package-test';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

category('splitters', async () => {
  before(async () => {
  });

  after(async () => {
  });

  const _helm1 = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$';

  const _helm2 = 'PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$';

  const data: { [key: string]: [string, string[]] } = {
    fastaMulti: [
      'M[MeI]YKETLL[MeF]PKTDFPMRGGL[MeA]',
      ['M', 'MeI', 'Y', 'K', 'E', 'T', 'L', 'L', 'MeF', 'P',
        'K', 'T', 'D', 'F', 'P', 'M', 'R', 'G', 'G', 'L', 'MeA'],
    ],
    fastaFromHelm: [
      '[meI][Pip][dK][Thr_PO3H2][L-hArg(Et,Et)][D-Tyr_Et][Tyr_ab-dehydroMe][dV]EN[D-Orn][D-aThr][Phe_4Me]',
      ['meI', 'Pip', 'dK', 'Thr_PO3H2', 'L-hArg(Et,Et)', 'D-Tyr_Et', 'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn',
        'D-aThr', 'Phe_4Me'],
    ],

    helm1: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et',
        'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn', 'D-aThr', 'Phe_4Me'],
    ],
    helm2: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.D-Dap.dV.E.N.pnG.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me'],
    ],
    // HELM editor dialog returns HELM string with multichar monomer names in square brackets
    helm3: [
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK].[Thr_PO3H2].[Aca].[D-Tyr_Et].[D-Dap].[dV].E.N.[pnG].[Phe_4Me]}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me'],
    ],

    testHelm1: [
      'RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$',
      ['R(U)P', 'R(T)P', 'R(G)P', 'R(C)P', 'R(A)'],
    ],

    testHelm2: [
      'RNA1{P.R(U)P.R(T)}$$$$',
      ['P', 'R(U)P', 'R(T)'],
    ],
    testHelm3: [
      'RNA1{P.R(U).P.R(T)}$$$$',
      ['P', 'R(U)', 'P', 'R(T)'],
    ],
  };

  test('fastaMulti', async () => { await _testFastaSplitter(data.fastaMulti[0], data.fastaMulti[1]); });
  test('fastaFromHelm', async () => { await _testFastaSplitter(data.fastaFromHelm[0], data.fastaFromHelm[1]); });

  test('helm1', async () => { await _testHelmSplitter(data.helm1[0], data.helm1[1]); });
  test('helm2', async () => { await _testHelmSplitter(data.helm2[0], data.helm2[1]); });
  test('helm3-multichar', async () => { await _testHelmSplitter(data.helm3[0], data.helm3[1]); });

  // examples from Helm/tests/test.csv file
  test('testHelm1', async () => { await _testHelmSplitter(data.testHelm1[0], data.testHelm1[1]); });
  test('testHelm2', async () => { await _testHelmSplitter(data.testHelm2[0], data.testHelm2[1]); });
  test('testHelm3', async () => { await _testHelmSplitter(data.testHelm3[0], data.testHelm3[1]); });
});

category('splitters', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  test('splitToMonomers', async () => {
    const df: DG.DataFrame = await grok.dapi.files.readCsv('System:AppData/Bio/samples/MSA.csv');

    const seqCol = df.getCol('MSA');
    const semType = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;
    seqCol.setTag(bioTAGS.aligned, C.MSA);

    const newDf = await splitToMonomersUI(df, seqCol);
    expect(newDf.columns.names().includes('17'), true);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(newDf);

    // TODO: Check cell.renderer for columns of monomers
    const tv: DG.TableView = grok.shell.addTableView(newDf);
    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);
  });

  test('getHelmMonomers', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(
      `HELM,Activity
PEPTIDE1{hHis.N.T}$$$,5.30751
PEPTIDE1{hHis.Aca.Cys_SEt}$$$,5.72388
`);
    await grok.data.detectSemanticTypes(df);
    const expectedMonomerList = ['hHis', 'Aca', 'Cys_SEt', 'N', 'T'];

    const helmCol: DG.Column = df.getCol('HELM');
    const res = getHelmMonomers(helmCol);

    const missed = expectedMonomerList.filter((m) => !res.includes(m));
    const unexpected = res.filter((m) => !expectedMonomerList.includes(m));
    if (missed.length > 0 || unexpected.length) {
      const msgs = [];
      if (missed.length > 0)
        msgs.push(`Missed monomers ${JSON.stringify(missed)}.`);
      if (unexpected.length > 0)
        msgs.push(`Unexpected monomers ${JSON.stringify(unexpected)}.`);

      throw new Error(msgs.join(' '));
    }
  });

  // test('helmAsFasta', async () => {
  //   // The columns can't be empty for SeqHandler
  //   /* eslint-disable max-len */
  //   const srcSeq = '[meI][Pip][dK][Thr_PO3H2][L-hArg(Et,Et)][D-Tyr_Et][Tyr_ab-dehydroMe][dV]EN[D-Orn][D-aThr][Phe_4Me]';
  //   const tgtSeqA = ['meI', 'Pip', 'dK', 'Thr_PO3H2', 'L-hArg(Et,Et)', 'D-Tyr_Et', 'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn', 'D-aThr', 'Phe_4Me'];
  //   /* eslint-enable max-len */
  //   const resSeqA = splitterAsFasta(srcSeq);
  //   expectArray(resSeqA, tgtSeqA);
  // });
});

export async function _testFastaSplitter(src: string, tgt: string[]) {
  const resSS: ISeqSplitted = splitterAsFasta(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(resSS)} .`);
  expectArray(wu.count(0).take(resSS.length).map((p) => resSS.getOriginal(p)).toArray(), tgt);
}

export async function _testHelmSplitter(src: string, tgt: string[]) {
  const resSS: ISeqSplitted = splitterAsHelm(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(resSS)} .`);
  expectArray(wu.count(0).take(resSS.length).map((p) => resSS.getOriginal(p)).toArray(), tgt);
}
