import * as DG from 'datagrok-api/dg';

import {getSimilarityFromDistance} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {_chemGetSimilarities, chemGetFingerprints} from '../chem-searches';
import {Fingerprint} from './chem-common';
import {DistanceMatrix, DistanceMatrixService, dmLinearIndex} from '@datagrok-libraries/ml/src/distance-matrix';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';

export async function getSimilaritiesMarix(dim: number, smiles: DG.Column, dfSmiles: DG.DataFrame,
  colName: string, simArr: (DG.Column | null)[]): Promise<(DG.Column | null)[]> {
  const fingerprints = await chemGetFingerprints(dfSmiles.col(colName)!, Fingerprint.Morgan, false);
  for (let i = 0; i != dim - 1; ++i) {
    const fp = fingerprints.shift();
    simArr[i] = fp ?
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances',
        _chemGetSimilarities(smiles.get(i)!, fingerprints)) : null;
  }
  return simArr;
}

export function getSimilaritiesMarixFromDistances(dim: number, distances: Float32Array, simArr: DG.Column[])
  : DG.Column[] {
  const linearIdx = dmLinearIndex(dim);
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = new Float32Array(dim - i - 1).fill(0);
    for (let j = i + 1; j < dim; ++j)
      similarityArr[j - i - 1] = getSimilarityFromDistance(distances[linearIdx(i, j)]);
    simArr[i] = DG.Column.fromFloat32Array('similarity', similarityArr);
  }
  return simArr;
}

/** The full pairwise Tanimoto similarity matrix of the molecules (`Similarity Matrix...`): a `symbol` column of
 * `symbols`' values, then a column per molecule. */
export async function similarityMatrix(molecules: DG.Column, symbols: DG.Column, fingerprintType: Fingerprint):
  Promise<DG.DataFrame> {
  const pb = DG.TaskBarProgressIndicator.create('Computing similarity matrix...');
  const dmService = new DistanceMatrixService(true, true);
  try {
    const fingerprints = await chemGetFingerprints(molecules, fingerprintType, false);
    const n = fingerprints.length;

    const symbolValues: string[] = symbols.toList().map((v) => v == null ? '' : v.toString());

    // parallel pairwise Tanimoto distances on web workers; condensed upper-triangular form
    const distances = await dmService.calc(fingerprints, BitArrayMetricsNames.Tanimoto, false);
    const dm = new DistanceMatrix(distances, n);

    const result = DG.DataFrame.create(n);
    result.name = `${molecules.name} similarity matrix`;
    result.columns.add(DG.Column.fromStrings('symbol', symbolValues));

    for (let j = 0; j < n; j++) {
      const values = new Float32Array(n);
      for (let i = 0; i < n; i++)
        values[i] = i === j ? 1 : 1 - dm.get(i, j);
      const colName = result.columns.getUnusedName(symbolValues[j] || `mol_${j}`);
      const col = DG.Column.fromFloat32Array(colName, values);
      col.meta.format = '0.000';
      result.columns.add(col);
    }
    return result;
  } finally {
    pb.close();
  }
}
