import * as DG from 'datagrok-api/dg';

import {sequenceChemSimilarity} from '../../monomer-works/monomer-utils';
import {ISeqSplitted} from '../../utils/macromolecule/types';
import {splitAlignedSequences} from '../splitter';
import {getSplitter} from './utils';
import {TAGS as bioTAGS} from '../../utils/macromolecule';
import {UnitsHandler} from '../units-handler';

export enum SCORE {
  IDENTITY = 'identity',
  SIMILARITY = 'similarity',
}

/** Performs transformations and sequence scoring.
 * @param {DG.DataFrame} table Table to attach results to.
 * @param {DG.Column<string>} col Sequences column to score. Must have Macromolecule semantic type.
 * @param {string} ref Reference sequence to score against.
 * @param {SCORE} scoring Scoring method.
 * @returns {DG.Column<number>} Scores column. */
export async function calculateScores(
  table: DG.DataFrame, col: DG.Column<string>, ref: string, scoring: SCORE
): Promise<DG.Column<number>> {
  const splitSeqDf = splitAlignedSequences(col);
  const srcUh = UnitsHandler.getOrCreate(col);
  const refCol = srcUh.getNewColumnFromList('ref', [ref]);
  const refUh = UnitsHandler.getOrCreate(refCol);
  const refSplitted = refUh.splitted[0]; // ref is at 0

  const scoresCol = scoring === SCORE.IDENTITY ? calculateIdentity(refSplitted, splitSeqDf) :
    scoring === SCORE.SIMILARITY ? await calculateSimilarity(refSplitted, splitSeqDf) : null;
  if (scoresCol === null)
    throw new Error(`In bio library: Unkown sequence scoring method: ${scoring}`);
  scoresCol.name = table.columns.getUnusedName(scoresCol.name);
  table.columns.add(scoresCol);
  return scoresCol;
}

/** Calculates identity scores as fraction of matching monomers on the same position.
 * @param {ISeqSplitted} reference Splitted reference sequence.
 * @param {DG.DataFrame} positionsDf Table which only contains position columns with semantic type Monomer.
 * @returns {DG.Column<number>} Scores column. */
export function calculateIdentity(reference: ISeqSplitted, positionsDf: DG.DataFrame): DG.Column<number> {
  const numPositions = positionsDf.columns.length;
  const positionCols: Uint32Array[] = new Array(numPositions);
  const positionEmptyCategories: number[] = new Array(numPositions);
  const categoryIndexesTemplate: number[] = new Array(numPositions);

  for (let posIdx = 0; posIdx < numPositions; ++posIdx) {
    const posCol = positionsDf.columns.byIndex(posIdx);
    positionCols[posIdx] = posCol.getRawData() as Uint32Array;
    positionEmptyCategories[posIdx] = posCol.categories.indexOf('');
    categoryIndexesTemplate[posIdx] = posCol.categories.indexOf(reference.getOriginal(posIdx) ?? '');
  }

  const identityScoresCol = DG.Column.float('Identity', positionsDf.rowCount);
  const identityScoresData = identityScoresCol.getRawData();
  for (let rowIndex = 0; rowIndex < positionsDf.rowCount; ++rowIndex) {
    identityScoresData[rowIndex] = 0;
    for (let posIdx = 0; posIdx < reference.length; ++posIdx) {
      const categoryIndex = positionCols[posIdx][rowIndex];
      if (categoryIndex === categoryIndexesTemplate[posIdx])
        ++identityScoresData[rowIndex];
    }
    identityScoresData[rowIndex] /= reference.length;
  }

  return identityScoresCol;
}

/** Calculates similarity scores as sum of monomer fingerprint similarities on the same position.
 * @param {ISeqSplitted} reference Splitted reference sequence.
 * @param {DG.DataFrame} positionsDf Table which only contains position columns with semantic type Monomer.
 * @return {DG.Column<number>} Scores column. */
export async function calculateSimilarity(
  reference: ISeqSplitted, positionsDf: DG.DataFrame
): Promise<DG.Column<number>> {
  const monomerColumns = positionsDf.columns.toList() as DG.Column<string>[];
  const scoresCol = await sequenceChemSimilarity(monomerColumns, reference);
  return scoresCol;
}
