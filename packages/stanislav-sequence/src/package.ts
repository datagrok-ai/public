/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: complement
//input: string nucleotides
//output: string result
export function complement(nucleotide: string): string {
  // your code goes here
  return nucleotide;
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