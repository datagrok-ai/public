import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {fuzzyJoin} from '../package';
import * as DG from 'datagrok-api/dg';
import {NUCLEOTIDE_SEMTYPE} from '../nucleotide-utils';

category('fuzzy join function', (): void => {
  const createDataFrameContainingDnaNucleotideCol = (colName: string, nucleotides: string[]): DG.DataFrame => {
    const nucleotidesCol = DG.Column.fromList(DG.TYPE.STRING, colName, nucleotides);
    nucleotidesCol.semType = NUCLEOTIDE_SEMTYPE;
    return DG.DataFrame.fromColumns([nucleotidesCol]);
  };

  test('check the existence of the appended columns and the "Counts" column, its values', async (): Promise<void> => {
    const N = 3;
    const countsColumnName = 'Counts';
    const seqColName = 'Sequence';

    const df1 = createDataFrameContainingDnaNucleotideCol(seqColName, [
      'aaa ttt',
      'ctt taa',
      'tgt gta',
    ]);

    const df2 = createDataFrameContainingDnaNucleotideCol(seqColName, [
      'ccc aaa',
      'ggg ttt',
      'tgt gaa',
    ]);

    const exptectedCountsValues = [2, 1, 3, 1, 2, 2];

    // Computations performing
    const dfResult = fuzzyJoin(df1, df2, N);

    // check dfResult to contain all the columns that the df1 has
    const df1ColumnsPresent = df1.columns.toList()
      .every((c: DG.Column<unknown>): boolean => dfResult.columns.contains(c.name));
    expect(df1ColumnsPresent, true);

    // check dfResult to contain the Counts column
    expect(dfResult.columns.contains(countsColumnName), true);

    // check Sequence column value join
    const expectedJoinedSeq = [...df1.getCol(seqColName).toList(), ...df2.getCol(seqColName).toList()];
    const actualJoinedSeq = dfResult.getCol(seqColName).toList();
    expectArray(actualJoinedSeq, expectedJoinedSeq);

    // check Counts column result values
    const actualCountsValues = dfResult.getCol(countsColumnName).toList();
    expectArray(actualCountsValues, exptectedCountsValues);
  });
});
