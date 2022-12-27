import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';

import {ConverterFunc} from './types';
import {NOTATION, NotationConverter, UnitsHandler} from '@datagrok-libraries/bio';

// import {mmSemType} from '../const';
// import {importFasta} from '../package';

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

    helmLoneDeoxyribose = 'helmLoneDeoxyribose',
    helmLoneRibose = 'helmLoneRibose',
    helmLonePhosphorus = 'helmLonePhosphorus',
    fastaLoneDeoxyribose = 'fastaLoneDeoxyribose',
    fastaLoneRibose = 'fastaLoneRibose',
    fastaLonePhosphorus = 'fastaLonePhosphorus',
  }

  const _csvTxts: { [key: string]: string } = {
    fastaPt: `seq
FWPHEY
YNRQWYV
MKPSEYV
`,
    separatorPt: `seq
F-W-P-H-E-Y
Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V
`,
    helmPt: `seq
PEPTIDE1{F.W.P.H.E.Y}$$$
PEPTIDE1{Y.N.R.Q.W.Y.V}$$$
PEPTIDE1{M.K.P.S.E.Y.V}$$$
`,
    fastaDna: `seq
ACGTC
CAGTGT
TTCAAC
`,
    separatorDna: `seq
A/C/G/T/C
C/A/G/T/G/T
T/T/C/A/A/C
`,
    helmDna: `seq
DNA1{D(A)P.D(C)P.D(G)P.D(T)P.D(C)P}$$$
DNA1{D(C)P.D(A)P.D(G)P.D(T)P.D(G)P.D(T)P}$$$
DNA1{D(T)P.D(T)P.D(C)P.D(A)P.D(A)P.D(C)P}$$$
`,
    fastaRna: `seq
ACGUC
CAGUGU
UUCAAC
`,
    separatorRna: `seq
A*C*G*U*C
C*A*G*U*G*U
U*U*C*A*A*C
`,
    helmRna: `seq
RNA1{R(A)P.R(C)P.R(G)P.R(U)P.R(C)P}$$$
RNA1{R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P}$$$
RNA1{R(U)P.R(U)P.R(C)P.R(A)P.R(A)P.R(C)P}$$$
`,
    fastaGaps: `seq
FW-PH-EYY
FYNRQWYV-
FKP-Q-SEYV
`,
    separatorGaps: `seq
F/W//P/H//E/Y/Y
F/Y/N/R/Q/W/Y/V/
F/K/P//Q//S/E/Y/V
`,
    helmGaps: `seq
PEPTIDE1{F.W.*.P.H.*.E.Y.Y}$$$
PEPTIDE1{F.Y.N.R.Q.W.Y.V.*}$$$
PEPTIDE1{F.K.P.*.Q.*.S.E.Y.V}$$$
`,
    helmLoneDeoxyribose: `seq
DNA1{D(A).D(C).D(G).D(T).D(C)}$$$
DNA1{D(C).D(A).D(G).D(T).D(G).D(T)P}$$$
DNA1{D(T).D(T).D(C).D(A).D(A).D(C)P}$$$
`,
    helmLoneRibose: `seq
RNA1{R(A).R(C).R(G).R(U).R(C)}$$$
RNA1{R(C).R(A).R(G).R(U).R(G).R(U)P}$$$
RNA1{R(U).R(U).R(C).R(A).R(A).R(C)P}$$$
`,
    helmLonePhosphorus: `seq
RNA1{P.P.R(A)P.R(C)P.R(G)P.R(U)P.R(C)P}$$$
RNA1{P.P.R(C)P.R(A)P.P.R(G)P.R(U)P.R(G)P.R(U)P}$$$
RNA1{P.R(U)P.R(U)P.R(C)P.R(A)P.R(A)P.R(C)P.P.P}$$$
`,
  };

  const _csvDfs: { [key: string]: Promise<DG.DataFrame> } = {};

  /** Also detects semantic types
   * @param {string} key
   * @return {Promise<DG.DataFrame>}
   */
  function readCsv(key: string): Promise<DG.DataFrame> {
    if (!(key in _csvDfs)) {
      _csvDfs[key] = (async (): Promise<DG.DataFrame> => {
        const csv: string = _csvTxts[key];
        const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
        await grok.data.detectSemanticTypes(df);
        return df;
      })();
    }
    return _csvDfs[key];
  };

  function converter(tgtNotation: NOTATION, tgtSeparator: string | null = null): ConverterFunc {
    if (tgtNotation === NOTATION.SEPARATOR && !tgtSeparator)
      throw new Error(`Argument 'separator' is missed for notation '${tgtNotation.toString()}'.`);

    return function(srcCol: DG.Column): DG.Column {
      const converter = new NotationConverter(srcCol);
      const resCol = converter.convert(tgtNotation, tgtSeparator);
      expect(resCol.getTag('units'), tgtNotation);
      return resCol;
    };
  };

  async function _testConvert(srcKey: string, converter: ConverterFunc, tgtKey: string) {
    const srcDf: DG.DataFrame = await readCsv(srcKey);
    const srcCol: DG.Column = srcDf.getCol('seq');

    // conversion results
    const resCol: DG.Column = converter(srcCol);

    // The correct reference data to compare conversion results with.
    const tgtDf: DG.DataFrame = await readCsv(tgtKey);
    const tgtCol: DG.Column = tgtDf.getCol('seq');

    expectArray(resCol.toList(), tgtCol.toList());
    const uh: UnitsHandler = new UnitsHandler(resCol);
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
