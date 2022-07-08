import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebLogo, SplitterFunc} from '@datagrok-libraries/bio/src/viewers/web-logo';

category('splitter', () => {
  const helm1 = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$';
  const helm2 = 'PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$';

  test('helm1', async () => { await _testHelmSplitter(helm1); });
  test('helm2', async () => { await _testHelmSplitter(helm2); });
});

export async function _testHelmSplitter(txt: string) {
  // const splitter: SplitterFunc = WebLogo.getSplitterAsHelm();
  //
  // const mList: string[] = splitter(txt);
  // expect(mList.length, 12);
}

