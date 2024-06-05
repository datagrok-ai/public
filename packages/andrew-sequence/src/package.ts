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
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function complement(nucleotides: string): string {
  const toReplace: Record<string, string> = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
  };
  const replaced: string[] = [];
  for (const char of nucleotides) {
    const upperChar = char.toUpperCase();
    replaced.push(upperChar in toReplace ? toReplace[upperChar] : char);
  }

  return replaced.join('');
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(value: string): DG.Widget {
  return new DG.Widget(value ?
    ui.divText('value - ' + value) :
    ui.divText('value is empty'),
  );
}
