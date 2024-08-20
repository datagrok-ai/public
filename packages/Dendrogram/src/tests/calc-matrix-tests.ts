import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  category,
  test,
  expectArray,
} from '@datagrok-libraries/utils/src/test';
import {mapToFixed} from './utils/array-utils';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TreeHelper} from '../utils/tree-helper';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
const sequences = [
  'CATGT',
  'AATGA',
  'AAGAT',
  'CCAGT',
];
const sequences2 = [
  'TATAG',
  'ATATA',
  'TAGAG',
  'GATTT',
];
const hNumbers = [1, 3, 5, 2];
const hNumbers2 = [3, 5, 2, 2.4];

const DNADistances1 = [0.4, 0.6, 0.4, 0.4, 0.8, 0.6];
const DNADistances2Cols = [0.707, 0.791, 0.901, 0.901, 1.25, 1.25];

const DNAAndNumericsDistances = [0.707, 1.25, 0.559, 0.707, 1.031, 1.061];
const numericsDistances = [2, 4, 1, 2, 1, 3];
const numericsDistancesTwoCols = [0.833, 1.054, 0.417, 1.118, 1.031, 0.75];

function setMacromoleculeTags(col: DG.Column) {
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = NOTATION.FASTA;
  col.setTag(TAGS.aligned, ALIGNMENT.SEQ);
  col.setTag(TAGS.alphabet, ALPHABET.DNA);
}

category('CalculateDistances', () => {
  test('CalcDistanceDNA', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCols = [DG.Column.fromStrings('Sequence', sequences)];
    const df = DG.DataFrame.fromColumns(seqCols);
    setMacromoleculeTags(seqCols[0]);
    const matrix = await th.calcDistanceMatrix(df, seqCols.map((col) => col.name));
    expectArray(mapToFixed(matrix!.data), mapToFixed(DNADistances1));
  });

  test('CalcDistanceDNATwoCols', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCol1 = DG.Column.fromStrings('Sequence', sequences);
    const seqCol2 = DG.Column.fromStrings('Sequence2', sequences2);
    setMacromoleculeTags(seqCol1);
    setMacromoleculeTags(seqCol2);
    const df = DG.DataFrame.fromColumns([seqCol1, seqCol2]);
    const matrix = await th.calcDistanceMatrix(df, ['Sequence', 'Sequence2']);
    expectArray(mapToFixed(matrix!.data), mapToFixed(DNADistances2Cols));
  });

  test('CalcDistanceDNAAndNumeric', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCol = DG.Column.fromStrings('Sequence', sequences);
    setMacromoleculeTags(seqCol);
    const numCol = DG.Column.fromList('int', 'numbers', hNumbers);
    const df = DG.DataFrame.fromColumns([seqCol, numCol]);
    const matrix = await th.calcDistanceMatrix(df, ['Sequence', 'numbers']);
    expectArray(mapToFixed(matrix!.data), mapToFixed(DNAAndNumericsDistances));
  });

  test('CalcDistanceNumeric', async () => {
    const th = new TreeHelper();
    const numCol = DG.Column.fromList('int', 'numbers', hNumbers);
    const df = DG.DataFrame.fromColumns([numCol]);
    const matrix = await th.calcDistanceMatrix(df, ['numbers']);
    expectArray(mapToFixed(matrix!.data), mapToFixed(numericsDistances));
  });

  test('CalcDistanceNumericTwoCols', async () => {
    const th = new TreeHelper();
    const numCol1 = DG.Column.fromList('int', 'numbers', hNumbers);
    const numCol2 = DG.Column.fromList('int', 'numbers2', hNumbers2);
    const df = DG.DataFrame.fromColumns([numCol1, numCol2]);
    const matrix = await th.calcDistanceMatrix(df, ['numbers', 'numbers2']);
    expectArray(mapToFixed(matrix!.data), mapToFixed(numericsDistancesTwoCols));
  });
});
