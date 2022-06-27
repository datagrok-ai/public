import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import { helmToFasta, helmToPeptide, helmToRNA } from '../package';
import { _package } from '../package-test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


category('Helm', () => {

  test('helmToFasta', async () => {
    expect(await helmToFasta('RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$'), '>RNA1UTGCA');
    expect(await helmToFasta('RNA1{P.R(U).P.R(T)}$$$$'), '>RNA1UT');
    expect(await helmToFasta('PEPTIDE1{A.G}$$$$V2.0'), '>PEPTIDE1AG')
  });

  test('helmToRNA', async () => {
    expect(await helmToRNA('RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$'), 'UTGCA');
    expect(await helmToRNA('RNA1{P.R(U).P.R(T)}$$$$'), 'UT');
    expect(await helmToRNA('RNA1{R(U)P.R(T)P}|RNA2{P.R(A)P.R(A)}$RNA1,RNA2,2:pair-6:pair|RNA1,RNA2,5:pair-3:pair$$$'), 'UT AA');
  });

  test('helmToPeptide', async () => {
    expect(await helmToPeptide('PEPTIDE1{A.G}$$$$V2.0'), 'AG');
    expect(await helmToPeptide('PEPTIDE1{L.V.A}|PEPTIDE2{L.V.A}$$$$'), 'LVA LVA');
    expect(await helmToPeptide('PEPTIDE1{A.R.C.A.A.K.T.C.D.A}$PEPTIDE1,PEPTIDE1,8:R3-3:R3$$$'), 'ARCAAKTCDA');
  });

  test('detectHelm', async () => {
    const file = await _package.files.readAsText('test.csv')
    const df = DG.DataFrame.fromCsv(file);
    let col = df.columns.byName('HELM string');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, DG.SEMTYPE.HELM);
  });
});
