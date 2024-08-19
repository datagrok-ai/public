import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {ConverterFunc} from './types';


category('converters', () => {
  enum Samples {
    fastaPt = 'fastaPt',
    separatorPt = 'separatorPt',
    helmPt = 'helmPt',

    fastaDna = 'fastaDna',
    separatorDna = 'separatorDna',
    helmDna = 'helmDna',

    fastaRna = 'fastaRna',
    separatorRna = 'separatorRna',
    helmRna = 'helmRna',

    fastaGaps = 'fastaGaps',
    separatorGaps = 'separatorGaps',
    helmGaps = 'helmGaps',

    fastaUn = 'fastaUn',
    separatorUn = 'separatorUn',
    helmUn = 'helmUn',

    helmLoneDeoxyribose = 'helmLoneDeoxyribose',
    helmLoneRibose = 'helmLoneRibose',
    helmLonePhosphorus = 'helmLonePhosphorus',
    fastaLoneDeoxyribose = 'fastaLoneDeoxyribose',
    fastaLoneRibose = 'fastaLoneRibose',
    fastaLonePhosphorus = 'fastaLonePhosphorus',
  }

  const _csvTxts: { [key: string]: string } = {
    [Samples.fastaPt]: `seq
FWPHEYFWPHEY
YNRQWYVYNRQWYV
MKPSEYVMKPSEYV`,
    [Samples.separatorPt]: `seq
F-W-P-H-E-Y-F-W-P-H-E-Y
Y-N-R-Q-W-Y-V-Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V-M-K-P-S-E-Y-V`,
    [Samples.helmPt]: `seq
PEPTIDE1{F.W.P.H.E.Y.F.W.P.H.E.Y}$$$$
PEPTIDE1{Y.N.R.Q.W.Y.V.Y.N.R.Q.W.Y.V}$$$$
PEPTIDE1{M.K.P.S.E.Y.V.M.K.P.S.E.Y.V}$$$$`,
    [Samples.fastaDna]: `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`,
    [Samples.separatorDna]: `seq
A/C/G/T/C/A/C/G/T/C
C/A/G/T/G/T/C/A/G/T/G/T
T/T/C/A/A/C/T/T/C/A/A/C`,
    [Samples.helmDna]: `seq
RNA1{d(A)p.d(C)p.d(G)p.d(T)p.d(C)p.d(A)p.d(C)p.d(G)p.d(T)p.d(C)p}$$$$
RNA1{d(C)p.d(A)p.d(G)p.d(T)p.d(G)p.d(T)p.d(C)p.d(A)p.d(G)p.d(T)p.d(G)p.d(T)p}$$$$
RNA1{d(T)p.d(T)p.d(C)p.d(A)p.d(A)p.d(C)p.d(T)p.d(T)p.d(C)p.d(A)p.d(A)p.d(C)p}$$$$`,
    [Samples.fastaRna]: `seq
ACGUCACGUC
CAGUGUCAGUGU
UUCAACUUCAAC`,
    [Samples.separatorRna]: `seq
A*C*G*U*C*A*C*G*U*C
C*A*G*U*G*U*C*A*G*U*G*U
U*U*C*A*A*C*U*U*C*A*A*C`,
    [Samples.helmRna]: `seq
RNA1{r(A)p.r(C)p.r(G)p.r(U)p.r(C)p.r(A)p.r(C)p.r(G)p.r(U)p.r(C)p}$$$$
RNA1{r(C)p.r(A)p.r(G)p.r(U)p.r(G)p.r(U)p.r(C)p.r(A)p.r(G)p.r(U)p.r(G)p.r(U)p}$$$$
RNA1{r(U)p.r(U)p.r(C)p.r(A)p.r(A)p.r(C)p.r(U)p.r(U)p.r(C)p.r(A)p.r(A)p.r(C)p}$$$$`,
    [Samples.fastaGaps]: `seq
FW-PH-EYYFW-PH-EYY
FYNRQWYV-FYNRQWYV-
FKP-Q-SEYVFKP-Q-SEYV`,
    [Samples.separatorGaps]: `seq
F/W//P/H//E/Y/Y/F/W//P/H//E/Y/Y
F/Y/N/R/Q/W/Y/V//F/Y/N/R/Q/W/Y/V/
F/K/P//Q//S/E/Y/V/F/K/P//Q//S/E/Y/V`,
    [Samples.helmGaps]: `seq
PEPTIDE1{F.W.*.P.H.*.E.Y.Y.F.W.*.P.H.*.E.Y.Y}$$$$
PEPTIDE1{F.Y.N.R.Q.W.Y.V.*.F.Y.N.R.Q.W.Y.V.*}$$$$
PEPTIDE1{F.K.P.*.Q.*.S.E.Y.V.F.K.P.*.Q.*.S.E.Y.V}$$$$`,

    [Samples.fastaUn]: `seq
[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]D[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]D
[meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2]
[Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca]`,
    [Samples.separatorUn]: `seq
meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-D-meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-D
meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2
Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca`,
    [Samples.helmUn]: `seq
PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D.meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D}$$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2}$$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca}$$$$`,
    [Samples.helmLoneDeoxyribose]: `seq
RNA1{d(A).d(C).d(G).d(T).d(C).d(A).d(C).d(G).d(T).d(C)}$$$$
RNA1{d(C).d(A).d(G).d(T).d(G).d(T)p.d(C).d(A).d(G).d(T).d(G).d(T)p}$$$$
RNA1{d(T).d(T).d(C).d(A).d(A).d(C)p.d(T).d(T).d(C).d(A).d(A).d(C)p}$$$$`,
    [Samples.helmLoneRibose]: `seq
RNA1{r(A).r(C).r(G).r(U).r(C).r(A).r(C).r(G).r(U).r(C)}$$$$
RNA1{r(C).r(A).r(G).r(U).r(G).r(U)p.r(C).r(A).r(G).r(U).r(G).r(U)p}$$$$
RNA1{r(U).r(U).r(C).r(A).r(A).r(C)p.r(U).r(U).r(C).r(A).r(A).r(C)p}$$$$`,
    [Samples.helmLonePhosphorus]: `seq
RNA1{p.p.r(A)p.r(C)p.r(G)p.r(U)p.r(C)p.r(A)p.r(C)p.r(G)p.r(U)p.r(C)p}$$$$
RNA1{p.p.r(C)p.r(A)p.p.r(G)p.r(U)p.r(G)p.r(U)p.r(C)p.r(A)p.p.r(G)p.r(U)p.r(G)p.r(U)p}$$$$
RNA1{p.r(U)p.r(U)p.r(C)p.r(A)p.r(A)p.r(C)p.r(U)p.r(U)p.r(C)p.r(A)p.r(A)p.r(C)p.p.p}$$$$`,
  };

  const bioTagsSet = new Set<string>(Object.values(bioTAGS));

  /** Also detects semantic types
   * @param {string} key
   * @return {Promise<DG.DataFrame>}
   */
  async function readCsv(key: string): Promise<DG.DataFrame> {
    // Always recreate test data frame from CSV for reproducible detector behavior in tests.
    const csv: string = _csvTxts[key];
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    await grok.data.detectSemanticTypes(df);
    return df;
  }

  function converter(tgtNotation: NOTATION, tgtSeparator?: string): ConverterFunc {
    if (tgtNotation === NOTATION.SEPARATOR && !tgtSeparator)
      throw new Error(`Argument 'separator' is mandatory for target notation '${tgtNotation.toString()}'.`);

    return function(srcCol: DG.Column): DG.Column {
      const converterSh = SeqHandler.forColumn(srcCol);
      const resCol = converterSh.convert(tgtNotation, tgtSeparator);
      expect(resCol.meta.units, tgtNotation);
      return resCol;
    };
  }

  async function _testConvert(srcKey: Samples, colConverter: ConverterFunc, tgtKey: Samples) {
    const srcDf: DG.DataFrame = await readCsv(srcKey);
    const srcCol: DG.Column = srcDf.getCol('seq');

    // conversion results
    const resCol: DG.Column = colConverter(srcCol);

    // The correct reference data to compare conversion results with.
    const tgtDf: DG.DataFrame = await readCsv(tgtKey);
    const tgtCol: DG.Column = tgtDf.getCol('seq');

    expectArray(resCol.toList(), tgtCol.toList());
    const srcSh: SeqHandler = SeqHandler.forColumn(srcCol);
    const resSh: SeqHandler = SeqHandler.forColumn(resCol);
    for (const [tagName, tgtTagValue] of Object.entries(tgtCol.tags)) {
      if (
        !bioTagsSet.has(tagName) ||
        (srcSh.notation === NOTATION.HELM && [bioTAGS.alphabet, bioTAGS.alphabetIsMultichar].includes(tagName as bioTAGS)) ||
        (resSh.notation === NOTATION.HELM && [bioTAGS.alphabet, bioTAGS.alphabetIsMultichar].includes(tagName as bioTAGS))
      ) continue;

      const resTagValue = resCol.getTag(tagName);
      expect(resTagValue, tgtTagValue,
        `Tag '${tagName}' expected value '${tgtTagValue}' is not equal to actual '${resTagValue}'.`);
    }
  }

  // FASTA tests
  // fasta -> separator
  test('testFastaPtToSeparator', async () => {
    await _testConvert(Samples.fastaPt, converter(NOTATION.SEPARATOR, '-'), Samples.separatorPt);
  });
  test('testFastaDnaToSeparator', async () => {
    await _testConvert(Samples.fastaDna, converter(NOTATION.SEPARATOR, '/'), Samples.separatorDna);
  });
  test('testFastaRnaToSeparator', async () => {
    await _testConvert(Samples.fastaRna, converter(NOTATION.SEPARATOR, '*'), Samples.separatorRna);
  });
  test('testFastaGapsToSeparator', async () => {
    await _testConvert(Samples.fastaGaps, converter(NOTATION.SEPARATOR, '/'), Samples.separatorGaps);
  });
  test('testFastaUnToSeparator', async () => {
    await _testConvert(Samples.fastaUn, converter(NOTATION.SEPARATOR, '-'), Samples.separatorUn);
  });

  // fasta -> helm
  test('testFastaPtToHelm', async () => {
    await _testConvert(Samples.fastaPt, converter(NOTATION.HELM), Samples.helmPt);
  });
  test('testFastaDnaToHelm', async () => {
    await _testConvert(Samples.fastaDna, converter(NOTATION.HELM), Samples.helmDna);
  });
  test('testFastaRnaToHelm', async () => {
    await _testConvert(Samples.fastaRna, converter(NOTATION.HELM), Samples.helmRna);
  });
  test('testFastaGapsToHelm', async () => {
    await _testConvert(Samples.fastaGaps, converter(NOTATION.HELM), Samples.helmGaps);
  });
  // TODO: testFastaUnToHelm
  // test('testFastaUnToHelm', async () => {
  //   await _testConvert(Samples.fastaUn, converter(NOTATION.HELM), Samples.helmUn);
  // });


  // SEPARATOR tests
  // separator -> fasta
  test('testSeparatorPtToFasta', async () => {
    await _testConvert(Samples.separatorPt, converter(NOTATION.FASTA), Samples.fastaPt);
  });
  test('testSeparatorDnaToFasta', async () => {
    await _testConvert(Samples.separatorDna, converter(NOTATION.FASTA), Samples.fastaDna);
  });
  test('testSeparatorRnaToFasta', async () => {
    await _testConvert(Samples.separatorRna, converter(NOTATION.FASTA), Samples.fastaRna);
  });
  test('testSeparatorGapsToFasta', async () => {
    await _testConvert(Samples.separatorGaps, converter(NOTATION.FASTA), Samples.fastaGaps);
  });
  test('testSeparatorUnToFasta', async () => {
    await _testConvert(Samples.separatorUn, converter(NOTATION.FASTA), Samples.fastaUn);
  });

  // separator -> helm
  test('testSeparatorPtToHelm', async () => {
    await _testConvert(Samples.separatorPt, converter(NOTATION.HELM), Samples.helmPt);
  });
  test('testSeparatorDnaToHelm', async () => {
    await _testConvert(Samples.separatorDna, converter(NOTATION.HELM), Samples.helmDna);
  });
  test('testSeparatorRnaToHelm', async () => {
    await _testConvert(Samples.separatorRna, converter(NOTATION.HELM), Samples.helmRna);
  });
  test('testSeparatorGapsToHelm', async () => {
    await _testConvert(Samples.separatorGaps, converter(NOTATION.HELM), Samples.helmGaps);
  });
  // TODO: testSeparatorUnToHelm
  // test('testSeparatorUnToHelm', async () => {
  //   await _testConvert(Samples.separatorUn, converter(NOTATION.HELM), Samples.helmUn);
  // });


  // HELM tests
  // helm -> fasta
  test('testHelmDnaToFasta', async () => {
    await _testConvert(Samples.helmDna, converter(NOTATION.FASTA), Samples.fastaDna);
  });
  test('testHelmRnaToFasta', async () => {
    await _testConvert(Samples.helmRna, converter(NOTATION.FASTA), Samples.fastaRna);
  });
  test('testHelmPtToFasta', async () => {
    await _testConvert(Samples.helmPt, converter(NOTATION.FASTA), Samples.fastaPt);
  });
  test('testHelmUnToFasta', async () => {
    await _testConvert(Samples.helmUn, converter(NOTATION.FASTA), Samples.fastaUn);
  });

  // helm -> separator
  test('testHelmDnaToSeparator', async () => {
    await _testConvert(Samples.helmDna, converter(NOTATION.SEPARATOR, '/'), Samples.separatorDna);
  });
  test('testHelmRnaToSeparator', async () => {
    await _testConvert(Samples.helmRna, converter(NOTATION.SEPARATOR, '*'), Samples.separatorRna);
  });
  test('testHelmPtToSeparator', async () => {
    await _testConvert(Samples.helmPt, converter(NOTATION.SEPARATOR, '-'), Samples.separatorPt);
  });
  test('testHelmUnToSeparator', async () => {
    await _testConvert(Samples.helmUn, converter(NOTATION.SEPARATOR, '-'), Samples.separatorUn);
  });

  // helm miscellaneous
  test('testHelmLoneRibose', async () => {
    await _testConvert(Samples.helmLoneRibose, converter(NOTATION.FASTA), Samples.fastaRna);
  });
  test('testHelmLoneDeoxyribose', async () => {
    await _testConvert(Samples.helmLoneDeoxyribose, converter(NOTATION.SEPARATOR, '/'), Samples.separatorDna);
  });
  test('testHelmLonePhosphorus', async () => {
    await _testConvert(Samples.helmLonePhosphorus, converter(NOTATION.FASTA), Samples.fastaRna);
  });
});
