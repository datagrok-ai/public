/**
 * The time limits of Bio's macromolecule assembly on the Crux workers (crux-service.ts `biopolymerChunkMs`): a chunk's
 * limit must cover the depictor's non-linear layout, whose time per monomer grows with the row (about 1-1.5 ms a
 * monomer at 250-400 monomers). A chunk past its limit restarts the pool and fails the whole column (reported:
 * non-linear To Atomic Level of 2,000 FASTA rows of 360-400 amino acids, 125 rows a chunk on 16 workers, threw "a Crux
 * job took longer than 78024 ms": 1 ms a monomer).
 */
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {AsmRole, AsmSequence, packBatch} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {biopolymerChunkMs} from '../crux/crux-service';

/** A peptide of `n` alanines. */
function chain(n: number): AsmSequence {
  return {monomers: Array.from({length: n}, (_, position) => ({symbol: 'A', polymerType: 'PEPTIDE' as const,
    role: AsmRole.AA, position})), graphInfo: {connections: [], disjointSeqStarts: n > 0 ? [0] : []}};
}

const allRows = (n: number): number[] => Array.from({length: n}, (_, i) => i);

category('crux: biopolymer: limits', () => {
  test('a chunk\'s limit covers the depictor on long rows', async () => {
    // the reported chunk: 125 rows of ~380 monomers; at least 2 ms a monomer (1.5 measured, and a loaded machine)
    const long = packBatch(Array.from({length: 125}, () => chain(380)));
    const longMs = biopolymerChunkMs(long, allRows(125));
    expect(longMs >= 2 * 125 * 380, true, `the limit of 125 rows of 380 monomers: ${longMs} ms`);
    // short rows keep the limit they had: 5 ms a row and 1 ms a monomer
    const short = packBatch(Array.from({length: 100}, () => chain(20)));
    expect(biopolymerChunkMs(short, allRows(100)), 100 * (5 + 20), 'the limit of 100 rows of 20 monomers');
    // a chunk of some rows of a batch counts those rows only
    const mixed = packBatch([chain(10), chain(300), chain(10)]);
    expect(biopolymerChunkMs(mixed, [0, 2]), 2 * (5 + 10), 'the limit of two of the rows');
  });
});
