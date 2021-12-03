/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
// eslint-disable-next-line no-unused-vars
import * as grok from 'datagrok-api/grok';
// eslint-disable-next-line no-unused-vars
import * as ui from 'datagrok-api/ui';
// eslint-disable-next-line no-unused-vars
import * as DG from 'datagrok-api/dg';

import biomsa from 'biomsa';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';

const BLOSUM62 = {
  matrix: [
    //    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    [4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2, 0], // A
    [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, 0], // C
    [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3, 0], // D
    [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2, 0], // E
    [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3, 0], // F
    [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3, 0], // G
    [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2, 0], // H
    [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1, 0], // I
    [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2, 0], // K
    [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1, 0], // L
    [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1, 0], // M
    [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2, 0], // N
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3, 0], // P
    [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1, 0], // Q
    [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2, 0], // R
    [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2, 0], // S
    [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2, 0], // T
    [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1, 0], // V
    [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2, 0], // W
    [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7, 0], // Y
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // *
  ],
  gapOP: -11,
  gapEP: -1,
};

export async function doMSA(col: DG.Column, gapopen: number, gapextend: number): Promise<DG.Column> {
  const sequences = col.toList().map((v: string, _) => AlignedSequenceEncoder.clean(v).replace(/\-/g, ''));
  const aligned = await biomsa.align(sequences, {
    matrix: BLOSUM62.matrix,
    gapopen: gapopen,
    gapextend: gapextend,
    method: 'complete',
    type: 'amino',
  });
  grok.shell.info([col.getRawData()[0].toString().length, aligned[0].length]);
  return DG.Column.fromStrings('aligned', aligned);
}
