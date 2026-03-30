import { category, expect, test } from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import { countSubsequenceTableAugment } from '../package';

category('CountSubsequenceTableAugment', () => {
  test('adds count column', async () => {
    const df = DG.DataFrame.fromCsv(
      `sequence,id
fasta: ctaccagaga,1
fasta: attaaaggtt,2
fasta: gttctctacc,3`
    );

    await countSubsequenceTableAugment(df, df.col('sequence')!, 'acc');

    expect(df.columns.contains('N(acc)'), true);
    expect(df.rowCount, 3);
    expect(df.get('N(acc)', 0), 1);
    expect(df.get('N(acc)', 1), 0);
    expect(df.get('N(acc)', 2), 1);
  });

  test('supports overlapping occurrences', async () => {
    const df = DG.DataFrame.fromCsv(
      `sequence,id
AAAA,1
ATA,2`
    );

    await countSubsequenceTableAugment(df, df.col('sequence')!, 'AA');

    expect(df.columns.contains('N(AA)'), true);
    expect(df.get('N(AA)', 0), 3);
    expect(df.get('N(AA)', 1), 0);
  });

  test('works with different subsequence names', async () => {
    const df = DG.DataFrame.fromCsv(
      `sequence,id
ATGATC,1
AAAA,2`
    );

    await countSubsequenceTableAugment(df, df.col('sequence')!, 'A');

    expect(df.columns.contains('N(A)'), true);
    expect(df.get('N(A)', 0), 2);
    expect(df.get('N(A)', 1), 4);
  });
});