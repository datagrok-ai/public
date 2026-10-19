import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {fuzzyJoin} from '../package';

function sequences(name: string, values: string[]): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(name, values)]);
  df.col(name)!.semType = 'dna_nucleotide';
  return df;
}

category('FuzzyJoin', () => {
  test('appends rows and counts matches', async () => {
    const df1 = sequences('seq', ['fasta: ATGC', 'fasta: GGG']);
    const df2 = sequences('dna', ['fasta: TGCA']);
    const result = fuzzyJoin(df1, df2, 2);

    expect(result.rowCount, 3);
    expect(result.columns.names().join(','), 'seq,Counts');
    // ATGC has AT, TG, GC: TG and GC occur once in TGCA. GGG has no matches.
    // TGCA has TG, GC, CA: TG and GC occur once in ATGC, none in GGG.
    expect(result.col('Counts')!.toList().join(','), '2,0,2');
  });
}, {clear: true});
