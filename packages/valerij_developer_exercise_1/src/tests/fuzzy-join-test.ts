import { category, expect, test } from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import { fuzzyJoin } from '../package';

category('fuzzyJoin', () => {
  test('appends dataframes and adds Counts column', async () => {
    const df1 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: GATTACA,1
fasta: ATTCGGA,2`
    );

    const df2 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: TTTAGGC,3
fasta: GATTACA,4`
    );

    df1.col('sequence')!.semType = 'dna_nucleotide';
    df2.col('sequence')!.semType = 'dna_nucleotide';

    const result = fuzzyJoin(df1, df2, 3);

    expect(result.rowCount, 4);
    expect(result.columns.contains('sequence'), true);
    expect(result.columns.contains('id'), true);
    expect(result.columns.contains('Counts'), true);
  });

  test('keeps unioned rows from both dataframes', async () => {
    const df1 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: AAAA,1
fasta: CCCC,2`
    );

    const df2 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: TTTT,3`
    );

    df1.col('sequence')!.semType = 'dna_nucleotide';
    df2.col('sequence')!.semType = 'dna_nucleotide';

    const result = fuzzyJoin(df1, df2, 2);

    expect(result.rowCount, 3);
    expect(result.get('id', 0), 1);
    expect(result.get('id', 1), 2);
    expect(result.get('id', 2), 3);
  });

  test('fills Counts with numeric values', async () => {
    const df1 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: AAAA,1`
    );

    const df2 = DG.DataFrame.fromCsv(
      `sequence,id
fasta: AAAA,2`
    );

    df1.col('sequence')!.semType = 'dna_nucleotide';
    df2.col('sequence')!.semType = 'dna_nucleotide';

    const result = fuzzyJoin(df1, df2, 2);
    const counts = result.col('Counts');

    expect(counts != null, true);
    expect(typeof result.get('Counts', 0), 'number');
    expect(typeof result.get('Counts', 1), 'number');
    expect(result.get('Counts', 0) > 0, true);
    expect(result.get('Counts', 1) > 0, true);
  });

  test('throws when dna_nucleotide column is missing', async () => {
    const df1 = DG.DataFrame.fromCsv(
      `text,id
AAAA,1`
    );

    const df2 = DG.DataFrame.fromCsv(
      `text,id
TTTT,2`
    );

    let errorMessage = '';

    try {
      fuzzyJoin(df1, df2, 2);
    } catch (e: any) {
      errorMessage = String(e.message ?? e);
    }

    expect(errorMessage.includes('dna_nucleotide'), true);
  });
});