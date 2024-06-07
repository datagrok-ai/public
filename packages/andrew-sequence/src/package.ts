/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as nu from './nucleotide-utils';

export const _package = new DG.Package();

//name: info
export function info(): void {
  grok.shell.info(_package.webRoot);
}

//name: complement
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function complement(nucleotides: string): string {
  return nu.complement(nucleotides);
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
  return await grok.functions.call(`${_package.name}:CountSubsequencePython`, {sequence, subsequence});
}

//name: callCountSubsequenceTableAugmentScript
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
export async function callCountSubsequenceTableAugmentScript(
  sequences: DG.DataFrame, columnName: DG.Column, subsequence: string,
): Promise<void> {
  const scriptName = `${_package.name}:CountSubsequencePythonDataframe`;
  const df = await grok.functions.call(scriptName, {sequences, columnName, subsequence});
  const countCol = df.columns.byIndex(0);
  countCol.name = `N(${subsequence})`;
  sequences.columns.insert(countCol);
}

//name: getOrders
//output: dataframe df
//input: string country = "USA"
export async function getOrders(country: string): Promise<DG.DataFrame> {
  return await grok.data.query(`${_package.name}:ordersByCountry`, {country});
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

//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: int N = 3
//output: dataframe result
export function fuzzyJoin(df1: DG.DataFrame, df2: DG.DataFrame, N: number): DG.DataFrame {
  const df = df1.append(df2);
  const countCol = df.columns.addNew('Counts', DG.TYPE.INT);
  const subsequencesCol = df.columns.bySemType(nu.NUCLEOTIDE_SEMTYPE);
  if (!subsequencesCol) return df;

  const df1Size = df1.rowCount;
  if (!(df1Size && df2.rowCount)) return df;

  const subsequenceLists = subsequencesCol.toList().map((seq: string): string[] => (
    nu.generateSubsequences(seq.replaceAll(/\s+/g, ''), N)
  ));

  for (let i = 0; i < subsequenceLists.length; ++i) {
    const isFirstfHalf = i < df1Size;
    const jEnd = isFirstfHalf ? subsequenceLists.length : df1Size;
    let matches = 0;
    for (let j = isFirstfHalf ? df1Size : 0; j < jEnd; ++j) {
      for (const seqEl of subsequenceLists[i])
        matches += Number(subsequenceLists[j].includes(seqEl));
    }
    countCol.set(i, matches);
  }

  return df;
}
