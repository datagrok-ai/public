import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {runPepsea} from '../utils/pepsea';
import {TestLogger} from './utils/test-logger';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {awaitContainerStart} from '../utils/docker';

category('PepSeA', () => {
  const testCsv = `HELM,MSA
"PEPTIDE1{F.L.R.G.W.[MeF].Y.S.N.N.C}$$$$","F.L.R.G.W.MeF.Y..S.N.N.C"
"PEPTIDE1{F.L.R.G.Y.[MeF].Y.W.N.C}$$$$","F.L.R.G.Y.MeF.Y.W...N.C"
"PEPTIDE1{F.G.Y.[MeF].Y.W.S.D.N.C}$$$$","F...G.Y.MeF.Y.W.S.D.N.C"
"PEPTIDE1{F.L.R.G.Y.[MeF].Y.W.S.N.D.C}$$$$","F.L.R.G.Y.MeF.Y.W.S.N.D.C"
"PEPTIDE1{F.V.R.G.Y.[MeF].Y.W.S.N.C}$$$$","F.V.R.G.Y.MeF.Y.W.S..N.C"
`;

  const pepseaStderrCsv: string = `HELM,MSA
"PEPTIDE1{F.L.Mis.G.W.[MeF].Y.S.N.N.C}$$$$","F.L.Mis.G.W.MeF.Y..S.N.N.C"
"PEPTIDE1{F.L.Mis.G.Y.[MeF].Y.W.N.C}$$$$","F.L.Mis.G.Y.MeF.Y...W.N.C"
"PEPTIDE1{F.G.Y.[MeF].Y.W.S.D.N.C}$$$$","F...G.Y.MeF.Y.W.S.D.N.C"
`;
  const pepseaStderrWarningList: string = 'Mis not found in Monomer Map\nMeF not found in Monomer Map\n';

  const pepseaErrorCsv: string = `HELM
"PEPTIDE1{[NH2].*.A.Q.T.T.Y.K.N.Y.R.R.N.L.L.*.[COOH]}$$$$"
"PEPTIDE1{[NH2].M.A.N.T.T.Y.K.N.Y.R.N.N.L.L.*.[COOH]}$$$$"
"PEPTIDE1{[NH2].*.A.N.T.T.Y.K.C.Y.R.R.N.L.L.*.[COOH]}$$$$"
"PEPTIDE1{[NH2].*.A.N.T.T.Y.K.F.Y.R.R.N.L.L.*.[COOH]}$$$$"
`;
  const pepseaErrorError: string = 'PepSeA error: The pair (*,M) couldn\'t be found in the substitution matrix';

  before(async () => {
    await awaitContainerStart();
  });

  test('Basic alignment', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const resMsaCol = await runPepsea(df.getCol('HELM'), 'msa(HELM)');
    const tgtMsaCol = df.getCol('MSA');
    for (let i = 0; i < resMsaCol!.length; ++i)
      expect(resMsaCol!.get(i) == tgtMsaCol.get(i), true);
  }, {timeout: 60000 /* docker */, stressTest: true});

  test('stderr', async () => {
    const logger = new TestLogger();
    const df = DG.DataFrame.fromCsv(pepseaStderrCsv);
    const resMsaCol = await runPepsea(df.getCol('HELM'), 'msa(HELM)',
      undefined, undefined, undefined, undefined, logger);
    const tgtMsaCol = df.getCol('MSA');
    expectArray(resMsaCol!.toList(), tgtMsaCol.toList());
    expect(logger.warningList[0].message, pepseaStderrWarningList);
  }, {timeout: 60000 /* docker */, stressTest: true});

  test('error', async () => {
    const logger = new TestLogger();
    try {
      const df = DG.DataFrame.fromCsv(pepseaErrorCsv);
      const _resMsaCol = await runPepsea(df.getCol('HELM'), 'msa(HELM)',
        undefined, undefined, undefined, undefined, logger);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      logger.error(errMsg, undefined, errStack);
    }
    expect(logger.errorList[0].message, pepseaErrorError);
  });
});
