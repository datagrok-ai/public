import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '@datagrok-libraries/bio/src/viewers/web-logo';

category('splitters', () => {
  const helm1 = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$';

  const helm2 = 'PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$';

  const data: { [key: string]: [string, string[]] } = {
    helm1: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et',
        'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn', 'D-aThr', 'Phe_4Me']
    ],
    helm2: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.D-Dap.dV.E.N.pnG.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me']
    ],
    // HELM editor dialog returns HELM string with multichar monomer names in square brackets
    helm3: [
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK].[Thr_PO3H2].[Aca].[D-Tyr_Et].[D-Dap].[dV].E.N.[pnG].[Phe_4Me]}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me']
    ],

    testHelm1: [
      'RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$',
      ['R(U)P', 'R(T)P', 'R(G)P', 'R(C)P', 'R(A)']
    ],

    testHelm2: [
      'RNA1{P.R(U)P.R(T)}$$$$',
      ['P', 'R(U)P', 'R(T)']
    ],
    testHelm3: [
      'RNA1{P.R(U).P.R(T)}$$$$',
      ['P', 'R(U)', 'P', 'R(T)']
    ],
  };

  test('helm1', async () => { await _testHelmSplitter(data.helm1[0], data.helm1[1]); });
  test('helm2', async () => { await _testHelmSplitter(data.helm2[0], data.helm2[1]); });
  test('helm3-multichar', async () => { await _testHelmSplitter(data.helm3[0], data.helm3[1]); });

  // examples from Helm/tests/test.csv file
  test('testHelm1', async () => { await _testHelmSplitter(data.testHelm1[0], data.testHelm1[1]); });
  test('testHelm2', async () => { await _testHelmSplitter(data.testHelm2[0], data.testHelm2[1]); });
  test('testHelm3', async () => { await _testHelmSplitter(data.testHelm3[0], data.testHelm3[1]); });
});

export async function _testHelmSplitter(src: string, tgt: string[]) {
  const res: string[] = WebLogo.splitterAsHelm(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(res)} .`);
  expectArray(res, tgt);
}

