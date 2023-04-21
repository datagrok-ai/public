import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  category,
  test,
  expectArray,
} from '@datagrok-libraries/utils/src/test';
import { mapToFixed } from './utils/array-utils';
import { ITreeHelper } from '@datagrok-libraries/bio/src/trees/tree-helper';
import { TreeHelper } from '../utils/tree-helper';
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

const DNADistances1 = [2, 3, 2, 2, 4, 3];
const DNADistances2Cols = [0.333, 0.5, 0.667, 0.667, 1.202, 1.118];
const DNAAndNumericsDistances = [0.333, 1.118, 0, 0.333, 1, 0.833];
const numericsDistances = [2, 4, 1, 2, 1, 3];
const numericsDistancesTwoCols = [0.745, 1.054, 0.333, 1.054, 1, 0.667];

category('CalculateDistances', () => {
  test('CalcDistanceDNA', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCols = [DG.Column.fromStrings('Sequence', sequences)];
    const df = DG.DataFrame.fromColumns(seqCols);
    seqCols[0].semType = DG.SEMTYPE.MACROMOLECULE;
    const matrix = await th.calcDistanceMatrix(df, seqCols.map(col => col.name));
    expectArray(mapToFixed(matrix!.data), mapToFixed(DNADistances1));
  });

  test('CalcDistanceDNATwoCols', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCol1 = DG.Column.fromStrings('Sequence', sequences);
    const seqCol2 = DG.Column.fromStrings('Sequence2', sequences2);
    seqCol1.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol2.semType = DG.SEMTYPE.MACROMOLECULE;
    const df = DG.DataFrame.fromColumns([seqCol1, seqCol2]);
    const matrix = await th.calcDistanceMatrix(df, ['Sequence', 'Sequence2']);
    expectArray(mapToFixed(matrix!.data), mapToFixed(DNADistances2Cols));
  });

  test('CalcDistanceDNAAndNumeric', async () => {
    const th: ITreeHelper = new TreeHelper();
    const seqCol = DG.Column.fromStrings('Sequence', sequences);
    const numCol = DG.Column.fromList('int', 'numbers', hNumbers);
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
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
