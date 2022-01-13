import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/bio/src/sequence-encoder';

// @ts-ignore
import initKalign from '../wasm/kalign';

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

export async function runKalign(col: DG.Column, webRootValue: string) : Promise<DG.Column> {
  // eslint-disable-next-line max-len
  // export EMCC_CFLAGS="-s ENVIRONMENT='web' -s MODULARIZE=1 -s 'EXTRA_EXPORTED_RUNTIME_METHODS=[\"FS\", \"callMain\"]' -s EXPORT_NAME='initKalign' -s FORCE_FILESYSTEM=1" && emmake make
  const _webPath = `${webRootValue}dist/kalign.wasm`;
  const sequences = col.toList().map((v: string, _) => AlignedSequenceEncoder.clean(v).replace(/\-/g, ''));
  const fasta = _stringsToFasta(sequences);

  const args = '-i input.fa -o output.fa'.split(' ');
  const kalign = await initKalign({
    noInitialRun: true,
    locateFile: () => _webPath,
  });
  kalign.FS.writeFile('input.fa', fasta, {encoding: 'utf8'});
  const ret = kalign.callMain(args);
  console.log(`main(argc, argv) returned ${ret}`);
  const buf = kalign.FS.readFile('output.fa', {encoding: 'utf8'}) as string;
  const aligned = _fastaToStrings(buf).slice(0, sequences.length);
  return DG.Column.fromStrings('aligned', aligned);
}
