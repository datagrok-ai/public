/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as nu from './nucleotide-utils';
import * as ena from './ena';
import * as fu from './file-units';
import * as su from './script-utils';

import NucleotideBoxCellRenderer from './nucleotide-cell-renderer';

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
    'value - ' + complement(value) :
    'value is empty',
  ));
}

//name: runCountSubsequenceScript
//input: string sequence {semType: dna_nucleotide}
//input: string subsequence
//input: string scriptType = "Python"
//output: int count
//test: runCountSubsequenceScript("A", "AAA", "Python") == 0
//test: runCountSubsequenceScript("aBbabaB", "aB", "JS") == 2
//test: runCountSubsequenceScript("ararar", "ararar", "Python") == 1
export async function runCountSubsequenceScript(
  sequence: string, subsequence: string, scriptType: su.ScriptType,
): Promise<number> {
  return await su.countSubsequences(_package.name, {sequence, subsequence}, scriptType);
}

//name: runCountSubsequenceForTableScript
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
//input: string scriptType = "Python"
export async function runCountSubsequenceForTableScript(
  sequences: DG.DataFrame, column: DG.Column, subsequence: string, scriptType: su.ScriptType,
): Promise<void> {
  const scriptPrams = {sequences, columnName: column, subsequence};
  const countCol = await su.getSubsequencesCountColumn(_package.name, scriptPrams, scriptType);
  countCol.name = `N(${subsequence})`;
  sequences.columns.insert(countCol);
}

//name: getOrders
//output: dataframe df
//input: string country = "USA"
export async function getOrders(country: string): Promise<DG.DataFrame> {
  return await su.selectOredersByCountry(_package.name, country);
}

//name: openTable
//input: string filepath
//input: string stategy = "dataFiles"
//output: dataframe df
export async function openTable(filepath: string, srategy: keyof fu.FileLoadStategies): Promise<DG.DataFrame> {
  const loader = fu.fileLoadStategies[srategy] ?? fu.fileLoadStategies.dataFiles;
  const df = await loader(filepath);
  grok.shell.addTableView(df);
  return df;
}

//name: Add Tables
export async function addTables(): Promise<void> {
  return await fu.addTablesFromPacakgeCsv(_package);
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

  // nu.conutMatchedBySeparatorArrayStrategy(subsequencesCol.toList(), N, df1Size)
  //   .forEach((m: number, i: number): void => countCol.set(i, m));

  nu.conutMatchedBySeparatorMapStrategy(subsequencesCol.toList(), N, df1Size)
    .forEach((m: number, i: number): void => countCol.set(i, m));

  return df;
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_nucleotide
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer(): NucleotideBoxCellRenderer {
  return new NucleotideBoxCellRenderer();
}

//name: ENA Sequence
//tags: panel, widgets
//input: string cellText {semType: EnaID}
//output: widget result
//condition: true
export async function enaSequence(cellText: string): Promise<DG.Widget> {
  const url = `${ena.ENA_API}/fasta/${cellText}`;
  const fasta = await grok.dapi.fetchProxy(url)
    .then((res: Response): Promise<string> => res.text());
  const enaSequence = ena.parseENASequenceFasta(fasta);

  return new DG.Widget(ui.box(
    ui.divV([
      ui.h1(enaSequence?.id ?? 'No ENA id'),
      ui.divText(enaSequence?.sequence ?? 'Sequence is empty'),
    ]),
  ));
}

//name: _fetchENASequence
//input: string query
//input: int limit = 3
//input: int offset = 0
//output: dataframe result
export async function _fetchENASequence(query: string, limit: number, offset: number): Promise<DG.DataFrame> {
  const {ids, sequences} = await ena.fetchSequences(query, limit, offset);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ids),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', sequences),
  ]);
  df.getCol('Sequence').semType = nu.NUCLEOTIDE_SEMTYPE;
  return df;
}

//name: formENADataTable
export function formENADataTable(): void {
  const previewDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', []),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', []),
  ]);
  previewDf.getCol('Sequence').semType = nu.NUCLEOTIDE_SEMTYPE;

  const fetchSequence = (params: {
    query: string,
    limit: number,
    offset: number,
  }) => _fetchENASequence(params.query, params.limit, params.offset);

  ena.showEnaSequenceFormDialog(previewDf, fetchSequence);
}
