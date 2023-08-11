
import {sequenceChemSimilarity} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import * as DG from 'datagrok-api/dg';

export function calculateIdentity(template: ISeqSplitted, splitSeqDf: DG.DataFrame): DG.Column<number> {
  const numPositions = splitSeqDf.columns.length;
  const positionCols: Uint32Array[] = new Array(numPositions);
  const positionEmptyCategories: number[] = new Array(numPositions);
  const categoryIndexesTemplate: number[] = new Array(numPositions);

  for (let posIdx = 0; posIdx < numPositions; ++posIdx) {
    const posCol = splitSeqDf.columns.byIndex(posIdx);
    positionCols[posIdx] = posCol.getRawData() as Uint32Array;
    positionEmptyCategories[posIdx] = posCol.categories.indexOf('');
    categoryIndexesTemplate[posIdx] = posCol.categories.indexOf(template[posIdx] ?? '');
  }

  const identityScoresCol = DG.Column.float('Identity', splitSeqDf.rowCount);
  const identityScoresData = identityScoresCol.getRawData();
  for (let rowIndex = 0; rowIndex < splitSeqDf.rowCount; ++rowIndex) {
    identityScoresData[rowIndex] = 0;
    for (let posIdx = 0; posIdx < template.length; ++posIdx) {
      const categoryIndex = positionCols[posIdx][rowIndex];
      if (categoryIndex === categoryIndexesTemplate[posIdx])
        ++identityScoresData[rowIndex];
    }
    identityScoresData[rowIndex] /= template.length;
  }

  return identityScoresCol;
}


export async function calculateSimilarity(template: ISeqSplitted, splitSeqDf: DG.DataFrame): Promise<DG.Column<number>> {
  const columns = splitSeqDf.columns.toList() as DG.Column<string>[];
  const scoresCol = await sequenceChemSimilarity(columns, template);
  return scoresCol;
}
