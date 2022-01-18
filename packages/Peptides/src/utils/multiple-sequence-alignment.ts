import * as DG from 'datagrok-api/dg';

//@ts-ignore
import Aioli from '@biowasm/aioli';

import {AlignedSequenceEncoder} from '@datagrok-libraries/bio/src/sequence-encoder';
//@ts-ignore
import {SEMTYPE} from '../semantics';

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
 * Converts aligned sequence to semantic type format.
 *
 * @param {string} seq Source sequence.
 * @return {string} Formatted sequence.
 */
function _castAligned(seq: string): string {
  let delimited = '';

  for (let i = 0; i < seq.length; ++i) {
    const char = seq[i];
    delimited += char == '-' ? char : `-${char}`;
  }
  return `NH2${delimited}-COOH`;
}

/**
 * Formats a batch of sequences to correspond the semantic type.
 *
 * @param {string[]} alignment List of aligned sequences.
 * @return {string[]} Formatted sequences.
 */
function _stringsToAligned(alignment: string[]): string[] {
  const nItems = alignment.length;
  const aligned = new Array<string>(nItems);

  for (let i = 0; i < nItems; ++i) {
    aligned[i] = _castAligned(alignment[i]);
  }
  return aligned;
}

/**
 * Runs Aioli environment with kalign tool.
 *
 * @param {DG.Column} col Column with sequences.
 * @param {boolean} isAligned Whether the column is aligned.
 * @return {Promise<DG.Column>} Aligned sequences.
 */
export async function runKalign(col: DG.Column, isAligned = false) : Promise<DG.Column> {
  let sequences = col.toList();

  if (isAligned) {
    sequences = sequences.map((v: string, _) => AlignedSequenceEncoder.clean(v).replace(/\-/g, ''));
  }

  const fasta = _stringsToFasta(sequences);

  const CLI = await new Aioli('kalign/3.3.1');
  await CLI.fs.writeFile('input.fa', fasta);
  const output = await CLI.exec(`kalign input.fa -f fasta -o result.fasta`);
  const buf = await CLI.cat('result.fasta');

  console.warn(output);

  const aligned = _fastaToStrings(buf).slice(0, sequences.length);
  const alignedCol = DG.Column.fromStrings(`(${col.name})msa`, _stringsToAligned(aligned));
  alignedCol.semType = SEMTYPE.ALIGNED;
  return alignedCol;
}
