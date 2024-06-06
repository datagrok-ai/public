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

//name: getOrders
//output: dataframe df
//input: string country = "USA"
export async function getOrders(country: string): Promise<DG.DataFrame> {
  return await grok.data.query('AndrewSequence:ordersByCountry', {country});
}

//input: string filepath
//output: dataframe df
export async function openTable1(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.getDemoTable(filepath);
  grok.shell.addTableView(df);
  return df;
}
//input: string filepath
//output: dataframe df
export async function openTable2(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.files.openTable(`System.DemoFiles/${filepath}`);
  grok.shell.addTableView(df);
  return df;
}

//input: string filepath
//output: dataframe df
export async function openTable3(filepath: string): Promise<DG.DataFrame> {
  const [df] = await (grok.functions.eval(`OpenServerFile("System:DemoFiles/${filepath}")`));
  grok.shell.addTableView(df);
  return df;
}

//name: Add Tables
export async function addTables(): Promise<void> {
  // Recursively list package files
  const files = await _package.files.list('', true);

  // Filter files by extension
  const csvFiles = files.filter((f) => f.extension === 'csv');

  // Load every table and add a view for it
  for (const file of csvFiles) {
    // Alternative ways to read a table are:
    // const df = await grok.data.loadTable(`${_package.webRoot}${file.name}`);
    
    // const df = await _package.files.readCsv(file.name);
    const df = await grok.data.files.openTable(`System:AppData/${_package.name}/${file.fileName}`);
    grok.shell.addTableView(df);
  }
}
