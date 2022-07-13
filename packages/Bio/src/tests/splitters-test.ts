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
  };

  test('helm1', async () => { await _testHelmSplitter(data.helm1[0], data.helm1[1]); });
  test('helm2', async () => { await _testHelmSplitter(data.helm2[0], data.helm2[1]); });
});

export async function _testHelmSplitter(src: string, tgt: string[]) {
  const res: string[] = WebLogo.splitterAsHelm(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(res)} .`);
  expectArray(res, tgt);
}

