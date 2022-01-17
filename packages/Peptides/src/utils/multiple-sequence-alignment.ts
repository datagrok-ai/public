import * as DG from 'datagrok-api/dg';

//@ts-ignore
import Aioli from '@biowasm/aioli';

/**
 * Converts array of sequences into simple fasta string.
 *
 * @param {string[]} sequences Input list of sequences.
 * @return {string} Fasta-formatted string.
 */
function _stringsToFasta(sequences: string[]): string {
  return sequences.reduce((a, v, i) => a + `>sample${i + 1}\n${v}\n`, '');
}

/**
 * Extracts array of sequences from simple fasta string.
 *
 * @param {string} fasta Fasta-formatted string.
 * @return {string[]} Output list of sequences.
 */
function _fastaToStrings(fasta: string): string[] {
  return fasta.replace(/>sample\d+(\r\n|\r|\n)/g, '').split('\n');
}

/**
 * Runs Aioli environment with kalign tool.
 *
 * @param {DG.Column} col Column with sequences.
 * @return {Promise<DG.Column>} Aligned sequences.
 */
export async function runKalign(col: DG.Column) : Promise<DG.Column> {
  const sequences = col.toList();
  const fasta = _stringsToFasta(sequences);

  const CLI = await new Aioli('kalign/3.3.1');
  await CLI.fs.writeFile('input.fa', fasta);
  const output = await CLI.exec(`kalign input.fa -f fasta -o result.fasta`);
  const buf = await CLI.cat('result.fasta');
  console.warn(output);

  const aligned = _fastaToStrings(buf).slice(0, sequences.length);
  return DG.Column.fromStrings(`${col.name}aligned`, aligned);
}
