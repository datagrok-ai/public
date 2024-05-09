/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { swapTwoWays } from './swap';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: complement
//input: string nucleotides
//output: string result
export async function complement(nucleotide: string): Promise<string> {
  let result = await swapTwoWays(nucleotide, 'A', 'T');
  return result;
}

//name: dna_nucleotide 
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function dna_nucleotide(nucleotide: string): string {
  // your code goes here
  return nucleotide;
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(nucleotide: any): any {
  // your code goes here
  return new DG.Widget(ui.divText(dna_nucleotide(nucleotide)));;
}

//name: executeFunction
//tags: panel, widgets
//input: string sequence 
//input: string subsequence    
//output: int result

export async function executeCountSubsequence(sequence: string, subsequence: string): Promise<number> {
  return await grok.functions.call('StanislavSequence:CountSubsequencePython', { sequence, subsequence });
}