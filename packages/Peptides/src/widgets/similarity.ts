
import * as DG from 'datagrok-api/dg';

export function calculateIdentity(template: string[], splitSeqDf: DG.DataFrame): DG.Column<number> {
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
