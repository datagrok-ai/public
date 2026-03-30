/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: complement
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function complement(nucleotides: string): string {
  const cleaned = nucleotides
    .toUpperCase()
    .replace(/^FASTA:\s*/, '')
    .trim();

  const map: Record<string, string> = {
    A: 'T',
    T: 'A',
    G: 'C',
    C: 'G',
  };

  return cleaned
    .split('')
    .map((ch) => map[ch] ?? ch)
    .join('');
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(nucleotides: string): DG.Widget {
  return new DG.Widget(
    ui.divV([
      ui.h2('Complement'),
      ui.divText(complement(nucleotides)),
    ])
  );
}

//name: CountSubsequencePythonPackageTS
//input: string sequence
//input: string subsequence
//output: int count
export async function countSubsequencePythonPackageTS(sequence: string, subsequence: string): Promise<number> {
  return await grok.functions.call(
    'valerij_developer_exercise_1:CountSubsequencePythonlocal',
    { sequence, subsequence }
  );
}

//name: CountSubsequenceTableAugment
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
export async function countSubsequenceTableAugment(sequences: DG.DataFrame, columnName: DG.Column, subsequence: string): Promise<void> {
   const df = await grok.functions.call('valerij_developer_exercise_1:CountSubsequencePythonDataframe', { sequences, columnName, subsequence });
   const countCol = df.columns.byIndex(0);
   countCol.name = `N(${subsequence})`;
   sequences.columns.insert(countCol);
}

//name: getOrders
//output: dataframe df
export async function getOrders() {
  return await grok.data.query('valerij_developer_exercise_1:ordersByCountry', { country: 'USA' });
}

//name: openTableViaDemo
//input: string filepath
//output: dataframe df
export async function openTableViaDemo(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.getDemoTable(filepath);
  grok.shell.addTableView(df);
  return df;
}

//name: openTableViaFiles
//input: string filepath
//output: dataframe df
export async function openTableViaFiles(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.files.openTable(`System:/${filepath}`);
  grok.shell.addTableView(df);
  return df;
}

//name: openTableViaServerFile
//input: string filepath
//output: dataframe df
export async function openTableViaServerFile(filepath: string): Promise<DG.DataFrame> {
  const df = (await grok.functions.eval(`OpenServerFile("System:DemoFiles/${filepath}")`))[0];
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
      const df = await _package.files.readCsv(file.name);
      grok.shell.addTableView(df);
      // Alternative ways to read a table are:
      // const df = await grok.data.loadTable(`${_package.webRoot}${file.path}`);
      // const df = await grok.data.files.openTable(`System:AppData/${_package.name}/${file.fileName}`);
   }
}