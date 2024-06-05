/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Nucleotide} from './types';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
const COMPLEMENT_TABLE: Record<Nucleotide, Nucleotide> = {
  'A': 'T',
  'T': 'A',
  'G': 'C',
  'C': 'G',
};

//name: complement
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function complement(nucleotides: string): string {
  return Array.from(nucleotides).map((char: string): string => {
    const upperChar = char.toUpperCase();
    const replaceable = upperChar in COMPLEMENT_TABLE;
    if (!replaceable) return char;
    const isUpperCase = upperChar === char;
    const replaced = COMPLEMENT_TABLE[upperChar as Nucleotide];
    return isUpperCase ? replaced : replaced.toLowerCase();
  }).join('');
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(value: string): DG.Widget {
  return new DG.Widget(ui.divText(value ?
    'value - ' + value :
    'value is empty',
  ));
}

//name: callCountSubsequencePythonScript
//input: string sequence {semType: dna_nucleotide}
//input: string subsequence
//output: int count
//test: callCountSubsequencePythonScript("A", "AAA") == 0
//test: callCountSubsequencePythonScript("aBbabaB", "aB") == 2
//test: callCountSubsequencePythonScript("ararar", "ararar") == 1
export async function callCountSubsequencePythonScript(sequence: string, subsequence: string): Promise<number> {
  return await grok.functions.call('AndrewSequence:CountSubsequencePython', {sequence, subsequence});
}

//name: callCountSubsequenceTableAugmentScript
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
export async function callCountSubsequenceTableAugmentScript(
  sequences: DG.DataFrame, columnName: DG.Column, subsequence: string,
): Promise<void> {
  const scriptName = 'AndrewSequence:CountSubsequencePythonDataframe';
  const df = await grok.functions.call(scriptName, {sequences, columnName, subsequence});
  const countCol = df.columns.byIndex(0);
  countCol.name = `N(${subsequence})`;
  sequences.columns.insert(countCol);
}
