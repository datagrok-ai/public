import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {ALIGNMENT, ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_testNeg, _testPos} from './detectors-tests';
import {DfReaderFunc} from './types';


category('detectors:weak-and-likely', () => {
  const enum csvTests {
    fastaDnaWeak1 = 'fastaDnaWeak1',
    fastaDnaWeak1LikelyName = 'fastaDnaWeak1LikelyName',

    fastaRnaWeak1 = 'fastaRnaWeak1',
    fastaRnaWeak1LikelyName = 'fastaRnaWeak1LikelyName',

    fastaPtWeak1 = 'fastaPtWeak1',
    fastaPtWeak1LikelyName = 'fastaPtWeak1LikelyName',

    fastaUn1 = 'fastaUn1',
    fastaUn1LikelyName = 'fastaUn1LikelyName',
    fastaUnMsa1LikelyName = 'fastaUnMsa1LikelyName',
  }

  const csvData: { [name: string]: string } = {
    [csvTests.fastaDnaWeak1]: `id,colName
1,TTTTT
2,TTTTT
3,TTTTT
4,TTTTT
`,
    [csvTests.fastaDnaWeak1LikelyName]: `id,seq
1,TTTTT
2,TTTTT
3,TTTTT
4,TTTTT
`,
    [csvTests.fastaRnaWeak1]: `id,colName
1,UUUUU
2,UUUUU
3,UUUUU
4,UUUUU
`,
    [csvTests.fastaRnaWeak1LikelyName]: `id,seq
1,UUUUU
2,UUUUU
3,UUUUU
4,UUUUU
`,
    [csvTests.fastaPtWeak1]: `id,colName
1,SLSLSPGK
2,SLSLSPGK
3,SLSLSPGK
4,SLSLSPGK
`,
    [csvTests.fastaPtWeak1LikelyName]: `id,seq
1,SLSLSPGK
2,SLSLSPGK
3,SLSLSPGK
4,SLSLSPGK
`,
    [csvTests.fastaUn1]: `id,colName
1,word
2,other
3,some
4,another
`,
    [csvTests.fastaUn1LikelyName]: `id,seq
1,word
2,other
3,some
4,another
`,
    [csvTests.fastaUnMsa1LikelyName]: `id,seq
1,word
2,male
3,bare
4,core
`,
  };

  const readCsv: (key: csvTests) => DfReaderFunc = (key: keyof typeof csvData) => {
    return async () => {
      // Always recreate test data frame from CSV for reproducible detector behavior in tests.
      const csv: string = csvData[key];
      const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
      await grok.data.detectSemanticTypes(df);
      return df;
    };
  };

  test(csvTests.fastaDnaWeak1, async () => {
    await _testNeg(readCsv(csvTests.fastaDnaWeak1), 'colName');
  });
  test(csvTests.fastaDnaWeak1LikelyName, async () => {
    await _testPos(readCsv(csvTests.fastaDnaWeak1LikelyName), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.DNA, 4, false);
  });

  test(csvTests.fastaRnaWeak1, async () => {
    await _testNeg(readCsv(csvTests.fastaRnaWeak1), 'colName');
  });
  test(csvTests.fastaRnaWeak1LikelyName, async () => {
    await _testPos(readCsv(csvTests.fastaRnaWeak1LikelyName), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.RNA, 4, false);
  });

  test(csvTests.fastaPtWeak1, async () => {
    await _testNeg(readCsv(csvTests.fastaPtWeak1), 'colName');
  });
  test(csvTests.fastaPtWeak1LikelyName, async () => {
    await _testPos(readCsv(csvTests.fastaPtWeak1LikelyName), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.PT, 20, false);
  });

  test(csvTests.fastaUn1, async () => {
    await _testNeg(readCsv(csvTests.fastaUn1), 'colName');
  });
  test(csvTests.fastaUn1LikelyName, async () => {
    await _testPos(readCsv(csvTests.fastaUn1LikelyName), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.UN, 11, false);
  });
  test(csvTests.fastaUnMsa1LikelyName, async () => {
    await _testPos(readCsv(csvTests.fastaUnMsa1LikelyName), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 10, false);
  });
});
