import * as DG from 'datagrok-api/dg';

import {calculateScores, SCORE} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

export async function calculateScoresWithEmptyValues(
  table: DG.DataFrame, macromolecule: DG.Column, reference: string, scoring: SCORE, seqHelper: ISeqHelper,
): Promise<DG.Column<number>> {
  const scores = await calculateScores(table, macromolecule, reference, scoring, seqHelper);
  for (let i = 0; i < scores.length; i++) {
    if (macromolecule.isNone(i))
      scores.set(i, null, false);
  }
  return scores;
}
