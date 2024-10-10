import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectObject, expectArray, before} from '@datagrok-libraries/utils/src/test';
import {
  NOTATION, getAlphabetSimilarity, monomerToShort, pickUpPalette, splitterAsFasta, splitterAsHelm,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {Nucleotides, NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {UnknownSeqPalette} from '@datagrok-libraries/bio/src/unknown';
import {getStatsForCol} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {GAP_SYMBOL} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

/** GAP_SYMBOL */
const g: string = GAP_SYMBOL;

category('bio', () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  const csvDfN1: string = `seq
ACGTCT
CAGTGT
TTCAAC`;

  /** 2 - is an error monomer
   * This sequence set should be classified as nucleotides sequences.
   * Small error, not similar to amino acids.
   */
  const csvDfN1e: string = `seq
ACGTAT
CAGTTG
TTCG2C`;

  /** Pure amino acids sequence */
  const csvDfAA1: string = `seq
FWPHEYV
YNRQWYV
MKPSEYV`;

  /** A - alanine, G - glycine, T -= threonine, C - cysteine, W - tryptophan
   * This sequence set should be detected as amino acids more than nucleotides.
   */
  const _csvDfAA2: string = `seq
AGTCAT
AGTCGC
AGTCATW`;

  /** This sequence set should be recognized as unknown. */
  const csvDfX: string = `seq
XZJ{}2
5Z4733
3Z6></
675687`;

  // anonymous functions specified in test() registering must return Promise<any>
  test('testGetStatsHelm1', async () => {
    const csv = `seq
PEPTIDE1{meI}$$$$`;
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const seqCol: DG.Column = df.getCol('seq')!;
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.HELM;
    const stats = getStatsForCol(seqCol, 1, splitterAsHelm);

    expectObject(stats.freq, {
      'meI': 1,
    });
    expect(stats.sameLength, true);
  });

  test('testGetStatsN1', async () => { await _testGetStats(csvDfN1); });
  test('testGetAlphabetSimilarity', async () => { await _testGetAlphabetSimilarity(); });

  test('testPickupPaletteN1', async () => { await _testPickupPaletteN1(csvDfN1); });
  test('testPickupPaletteN1e', async () => { await _testPickupPaletteN1e(csvDfN1e); });
  test('testPickupPaletteAA1', async () => { await _testPickupPaletteAA1(csvDfAA1); });
  test('testPickupPaletteX', async () => { await _testPickupPaletteX(csvDfX); });

  function _testGetStats(csvDfN1: string) {
    const dfN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1);
    const seqCol: DG.Column = dfN1.col('seq')!;
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    const stats = getStatsForCol(seqCol, 5, splitterAsFasta);

    expectObject(stats.freq, {
      'A': 4,
      'C': 5,
      'G': 3,
      'T': 6,
    });
    expect(stats.sameLength, true);
  }

  async function _testGetAlphabetSimilarity() {
    const freq: { [m: string]: number } = {
      'A': 2041,
      'C': 3015,
      'G': 3015,
      'T': 2048,
      [g]: 1000,
    };
    const alphabet: Set<string> = new Set(Object.keys(Nucleotides.Names));
    const res = getAlphabetSimilarity(freq, alphabet);

    expect(res > 0.6, true);
  }

  async function _testPickupPaletteN1(csvDfN1: string) {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1);
    const col: DG.Column = df.col('seq')!;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    const cp = pickUpPalette(col, seqHelper);

    expect(cp instanceof NucleotidesPalettes, true);
  }

  async function _testPickupPaletteN1e(csvDfN1e: string) {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1e);
    const col: DG.Column = df.col('seq')!;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    const cp = pickUpPalette(col, seqHelper);

    expect(cp instanceof NucleotidesPalettes, true);
  }

  async function _testPickupPaletteAA1(csvDfAA1: string) {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDfAA1);
    const col: DG.Column = df.col('seq')!;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    const cp = pickUpPalette(col, seqHelper);

    expect(cp instanceof AminoacidsPalettes, true);
  }

  async function _testPickupPaletteX(csvDfX: string) {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDfX);
    const col: DG.Column = df.col('seq')!;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.FASTA;
    const cp = pickUpPalette(col, seqHelper);

    expect(cp instanceof UnknownSeqPalette, true);
  }

  async function _testPickupPaletteAA2(dfAA2: DG.DataFrame) {
    const seqCol: DG.Column = dfAA2.col('seq')!;
    const cp = pickUpPalette(seqCol, seqHelper);

    expect(cp instanceof AminoacidsPalettes, true);
  }
});

category('WebLogo.monomerToShort', () => {
  test('longMonomerSingle', async () => {
    expect(monomerToShort('S', 5), 'S');
  });
  test('longMonomerShort', async () => {
    expect(monomerToShort('Short', 5), 'Short');
  });
  test('longMonomerLong56', async () => {
    expect(monomerToShort('Long56', 6), 'Long56');
  });
  test('longMonomerComplexFirstPartShort', async () => {
    expect(monomerToShort('Long-long', 5), 'Long…');
  });
  test('longMonomerComplexFirstPartLong56', async () => {
    expect(monomerToShort('Long56-long', 6), 'Long5…');
  });
  test('monomerToShort', async () => {
    const pairs = [
      ['AbC', 'AbC'],
      ['AbCd', 'Ab…'],
      ['ABc', 'ABc'],
      ['ABcd', 'AB…'],
      ['A_b', 'A_b'],
      ['A_bc', 'A…'],
      ['Ab_c', 'Ab…'],
      ['A1_b', 'A1…'],
      ['Abc_d', 'Ab…'],
      ['Abcd_e', 'Ab…'],
      ['A-b', 'A-b'],
      ['A-bc', 'A…'],
      ['Ab-c', 'Ab…'],
      ['A1-b', 'A1…'],
      ['Abc-d', 'Ab…'],
      ['Abcd-e', 'Ab…'],
      ['A', 'A'],
      ['Ab', 'Ab'],
      ['Abc', 'Abc'],
      ['Ab…', 'Ab…'],
      ['Abcd', 'Ab…'],
      ['Abcde', 'Ab…'],
    ];
    const src: string[] = pairs.map((p) => p[0]);
    const tgt: string[] = pairs.map((p) => p[1]);
    const res: string [] = src.map((m) => monomerToShort(m, 3));
    expectArray(res, tgt);
  });
});
