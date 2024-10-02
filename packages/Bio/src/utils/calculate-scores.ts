import * as DG from 'datagrok-api/dg';
import { calculateScores, SCORE } from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';

export async function calculateScoresWithEmptyValues(
    table: DG.DataFrame, macromolecule: DG.Column, reference: string, scoring: SCORE
): Promise<DG.Column<number>> {
    const scores = await calculateScores(table, macromolecule, reference, scoring);
    for (let i = 0; i < scores.length; i++) {
        if (macromolecule.isNone(i))
            scores.set(i, null, false);
    }
    return scores;
}