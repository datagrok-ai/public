/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
//import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';
import {transposeMatrix} from '@datagrok-libraries/utils/src/operations';
import {Matrix, Vector} from '@datagrok-libraries/utils/src/type_declarations';

const Spearman = require('spearman-rho');

export function calcPositions(col: DG.Column): Matrix {
  const sequences = col.toList().map((v, _) => AlignedSequenceEncoder.clean(v));
  const enc = new AlignedSequenceEncoder();
  const encSeqs = sequences.map((v) => Vector.from(enc.encode(v)));
  const positions = transposeMatrix(encSeqs);
  return positions;
}

export async function calcSpearmanRhoMatrix(positions: Matrix): Promise<Matrix> {
  const nItems = positions.length;
  const rho = new Array(nItems).fill(0).map((_) => new Float32Array(nItems).fill(0));

  for (let i = 0; i < nItems; ++i) {
    for (let j = i+1; j < nItems; ++j) {
      rho[i][j] = await(new Spearman(positions[i], positions[j])).calc();
      rho[j][i] = rho[i][j];
    }
  }
  return rho;
}
