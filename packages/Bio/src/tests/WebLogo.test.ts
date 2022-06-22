import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Nucleotides, NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {Aminoacids, AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnknownSeqPalette} from '@datagrok-libraries/bio/src/unknown';

category('WebLogo', () => {
  test('testGetAlphabetSimilarity', async () => { _testGetAlphabetSimilarity(); });

  test('testPickupPaletteN1', async () => { _testPickupPaletteN1(); });
  test('testPickupPaletteAA1', async () => { _testPickupPaletteAA1(); });

  // dfAA2 is too similar to nucleotides
  // test('testPickupPaletteAA2', async () => {
  //   _testPickupPaletteAA2();
  // });
});

const dfN1: DG.DataFrame = DG.DataFrame.fromCsv(
  `seq
ACGTCT
CAGTGT
TTCAAC
`);

/** 2 - is an error monomer
 * This sequence set should be classified as nucleotides sequences.
 * Small error, not similar to amino acids.
 */
const dfN1e: DG.DataFrame = DG.DataFrame.fromCsv(
  `seq
ACGTAT
CAGTTG
TTCG2C
`);

/** Pure amino acids sequence */
const dfAA1: DG.DataFrame = DG.DataFrame.fromCsv(
  `seq
FWPHEYV
YNRQWYV
MKPSEYV
`);

/** A - alanine, G - glycine, T -= threonine, C - cysteine, W - tryptophan
 * This sequence set should be detected as amino acids more than nucleotides.
 */
const dfAA2: DG.DataFrame = DG.DataFrame.fromCsv(
  `seq
AGTCAT
AGTCGC
AGTCATW
`);

/** This sequence set should be recognized as unknown. */
const dfX: DG.DataFrame = DG.DataFrame.fromCsv(
  `seq
XZJ{}2
5Z4733
3Z6></
675687
`);

export function _testGetAlphabetFreqs() {
  const seqCol: DG.Column = dfN1.col('seq')!;
  const mFreq = WebLogo.getAlphabetFreqs(seqCol);

  expectObject(mFreq, {
    'A': 4,
    'C': 5,
    'G': 3,
    'T': 6
  });
}

export function _testGetAlphabetSimilarity() {
  const freq: { [m: string]: number } = {
    'A': 2041,
    'C': 3015,
    'G': 3015,
    'T': 2048,
    '-': 1000
  };
  const alphabet: Set<string> = new Set(Object.keys(Nucleotides.Names));
  const res = WebLogo.getAlphabetSimilarity(freq, alphabet);

  expect(res > 0.6, true);
}

export function _testPickupPaletteN1() {
  const seqCol: DG.Column = dfN1.col('seq')!;
  const cp = WebLogo.pickUpPalette(seqCol);

  expect(cp instanceof NucleotidesPalettes, true);
}

export function _testPickupPaletteAA1() {
  const seqCol: DG.Column = dfAA1.col('seq')!;
  const cp = WebLogo.pickUpPalette(seqCol);

  expect(cp instanceof AminoacidsPalettes, true);
}

export function _testPickupPaletteAA2() {
  const seqCol: DG.Column = dfAA2.col('seq')!;
  const cp = WebLogo.pickUpPalette(seqCol);

  expect(cp instanceof AminoacidsPalettes, true);
}

export function _testPickupPaletteAll() {
  const seqColN1: DG.Column = dfN1.col('seq')!;
  const seqColAA1: DG.Column = dfAA1.col('seq')!;
  const seqColAA2: DG.Column = dfAA2.col('seq')!;
  const seqColX: DG.Column = dfX.col('seq')!;

  const cpN1: SeqPalette = WebLogo.pickUpPalette(seqColN1);
  const cpAA1: SeqPalette = WebLogo.pickUpPalette(seqColAA1);
  const cpAA2: SeqPalette = WebLogo.pickUpPalette(seqColAA2);
  const cpX: SeqPalette = WebLogo.pickUpPalette(seqColX);

  expect(cpN1 instanceof NucleotidesPalettes, true);
  expect(cpAA1 instanceof AminoacidsPalettes, true);
  expect(cpAA2 instanceof AminoacidsPalettes, true);
  expect(cpX instanceof UnknownSeqPalette, true);
}
