/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {runKalign} from '../utils/multiple-sequence-alignment';
import {multipleSequenceAlignmentUI} from '../utils/multiple-sequence-alignment-ui';

category('MSA', async () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  //table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  const fromCsv = `seq
FWRWYVKHP
YNRWYVKHP
MWRSWYCKHP`;
  const toCsv = `seq
FWR-WYVKHP
YNR-WYVKHP
MWRSWYCKHP`;

  const longFromCsv = `seq
FWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHP
YNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHP
MWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHP`;

  const longToCsv = `seq
FWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHP
YNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHP
MWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHP`;


  const helmFromCsv = `seq
  PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2}$$$$
  PEPTIDE1{meI.Aca.N.T.dE.Thr_PO3H2}$$$$
  PEPTIDE1{hHis.Aca.N.T.dE.Thr_PO3H2}$$$$`;

  const helmToCsv = `seq
  meI.hHis.Aca.N.T.dE.Thr_PO3H2
  .meI.Aca.N.T.dE.Thr_PO3H2
  .hHis.Aca.N.T.dE.Thr_PO3H2`;

  const longHelmFromCsv = `seq
  PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.dV.E.N.D-Orn.D-aThr.Phe_4Me.Thr_PO3H2}$$$$
  PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Tyr_PO3H2.D-Chg.dV.Phe_ab-dehydro.N.D-Orn.D-aThr.Phe_4Me}$$$$
  PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Tyr_PO3H2.D-Chg.dV.Thr_PO3H2.N.D-Orn.D-aThr.Phe_4Me}$$$$`;

  const longHelmToCsv = `seq
  meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.dV.E.N.D-Orn.D-aThr.Phe_4Me.Thr_PO3H2
  meI.hHis.Aca.Cys_SEt.T.dK..Tyr_PO3H2.D-Chg.dV.Phe_ab-dehydro.N.D-Orn.D-aThr.Phe_4Me.
  Lys_Boc.hHis.Aca.Cys_SEt.T.dK..Tyr_PO3H2.D-Chg.dV.Thr_PO3H2.N.D-Orn.D-aThr.Phe_4Me.`;

  const SeparatorFromCsv = `seq
  F-W-P-H-E-Y
  Y-N-R-Q-W-Y-V
  M-K-P-S-E-Y-V`;

  const SeparatorToCsv = `seq
  FWPHEY-
  YNRQWYV
  MKPSEYV`;

  const SeparatorLongFromCsv = `seq
  M-I-E-V-F-L-F-G-I-V-L-G-L-I-P-I-T-L-A-G-L-F-V-T-A-Y-L-Q-Y-R-R-G-D-Q-L-D-L
  M-M-E-L-V-L-K-T-I-I-G-P-I-V-V-G-V-V-L-R-I-V-D-K-W-L-N-K-D-K
  M-D-R-T-D-E-V-S-N-H-T-H-D-K-P-T-L-T-W-F-E-E-I-F-E-E-Y-H-S-P-F-H-N`;

  const SeparatorLongToCsv = `seq
  MIEV-FLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL
  MMEL-VLKTII-GPIVVGVVLRIVDKWLNKDK------
  MDRTDEVSNHTHDKPTLTWFEEIFEEYHSPFHN-----`;

  test('isCorrect', async () => {
    await _testMsaIsCorrect(fromCsv, toCsv);
  });

  test('isCorrectLong', async () => {
    await _testMsaIsCorrect(longFromCsv, longToCsv);
  });

  test('isCorrectHelm', async () => {
    await _testMSAOnColumn(helmFromCsv, helmToCsv, NOTATION.HELM, NOTATION.SEPARATOR, undefined, 'mafft');
  }, {timeout: 80000 /* docker */});

  test('isCorrectHelmLong', async () => {
    await _testMSAOnColumn(longHelmFromCsv, longHelmToCsv, NOTATION.HELM, NOTATION.SEPARATOR, undefined, 'mafft');
  }, {timeout: 80000 /* docker */});

  test('isCorrectSeparator', async () => {
    await _testMSAOnColumn(
      SeparatorFromCsv, SeparatorToCsv, NOTATION.SEPARATOR, NOTATION.FASTA, ALPHABET.PT,
    );
  });

  test('isCorrectSeparatorLong', async () => {
    await _testMSAOnColumn(
      SeparatorLongFromCsv, SeparatorLongToCsv, NOTATION.SEPARATOR, NOTATION.FASTA, ALPHABET.PT,
    );
  });

  async function _testMSAOnColumn(
    srcCsv: string, tgtCsv: string,
    srcNotation: NOTATION, tgtNotation: NOTATION, alphabet?: ALPHABET, pepseaMethod?: string,
  ): Promise<void> {
    const srcDf: DG.DataFrame = DG.DataFrame.fromCsv(srcCsv);
    await grok.data.detectSemanticTypes(srcDf);
    const tgtDf: DG.DataFrame = DG.DataFrame.fromCsv(tgtCsv);

    const tgtCol = tgtDf.getCol('seq')!;
    const srcCol: DG.Column = srcDf.getCol('seq')!;
    expect(srcCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(srcCol.meta.units, srcNotation);
    if (alphabet)
      expect(srcCol.getTag(bioTAGS.alphabet), alphabet);

    const msaSeqCol = await multipleSequenceAlignmentUI({col: srcCol, pepsea: {method: pepseaMethod}}, seqHelper);
    expect(msaSeqCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(msaSeqCol.meta.units, tgtNotation);
    expect(msaSeqCol.getTag(bioTAGS.aligned), ALIGNMENT.SEQ_MSA);
    if (alphabet)
      expect(msaSeqCol.getTag(bioTAGS.alphabet), alphabet);
    expectArray(msaSeqCol.toList(), tgtCol.toList());
  }
});

async function _testMsaIsCorrect(srcCsv: string, tgtCsv: string): Promise<void> {
  const srcDf: DG.DataFrame = DG.DataFrame.fromCsv(srcCsv);
  await grok.data.detectSemanticTypes(srcDf);
  const tgtDf: DG.DataFrame = DG.DataFrame.fromCsv(tgtCsv);

  const srcCol: DG.Column = srcDf.getCol('seq')!;
  expect(srcCol.semType, DG.SEMTYPE.MACROMOLECULE);

  const tgtCol: DG.Column = tgtDf.getCol('seq')!;
  const resCol: DG.Column = await runKalign(srcCol, true);
  expectArray(resCol.toList(), tgtCol.toList());
}


