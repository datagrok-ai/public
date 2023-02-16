import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {runPepsea} from '../utils/pepsea';

category('PepSeA', () => {
  const testCsv = `HELM,MSA
  "PEPTIDE1{F.L.R.G.W.[MeF].Y.S.N.N.C}$$$$","F.L.R.G.W.MeF.Y..S.N.N.C"
  "PEPTIDE1{F.L.R.G.Y.[MeF].Y.W.N.C}$$$$","F.L.R.G.Y.MeF.Y.W...N.C"
  "PEPTIDE1{F.G.Y.[MeF].Y.W.S.D.N.C}$$$$","F...G.Y.MeF.Y.W.S.D.N.C"
  "PEPTIDE1{F.L.R.G.Y.[MeF].Y.W.S.N.D.C}$$$$","F.L.R.G.Y.MeF.Y.W.S.N.D.C"
  "PEPTIDE1{F.V.R.G.Y.[MeF].Y.W.S.N.C}$$$$","F.V.R.G.Y.MeF.Y.W.S..N.C"`;

  test('Basic alignment', async () => {
    const table = DG.DataFrame.fromCsv(testCsv);
    const alignedCol = await runPepsea(table.getCol('HELM'), 'msa(HELM)');
    const alignedTestCol = table.getCol('MSA');
    for (let i = 0; i < alignedCol.length; ++i)
      expect(alignedCol.get(i) == alignedTestCol.get(i), true);
  });
});
