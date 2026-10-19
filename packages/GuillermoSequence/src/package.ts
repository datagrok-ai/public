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

const COMPLEMENTS: {[n: string]: string} = {A: 'T', T: 'A', G: 'C', C: 'G', a: 't', t: 'a', g: 'c', c: 'g'};
const FASTA_PREFIX = /^fasta:\s*/;

//name: complement
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function complement(nucleotides: string): string {
  const prefix = nucleotides.match(FASTA_PREFIX)?.[0] ?? '';
  return prefix + nucleotides.slice(prefix.length).replace(/[ATGCatgc]/g, (n) => COMPLEMENTS[n]);
}

function normalize(sequence: string): string {
  return sequence.replace(FASTA_PREFIX, '').replace(/\s/g, '').toUpperCase();
}

function countSubsequences(sequences: string[], n: number): Map<string, number> {
  const counts = new Map<string, number>();
  for (const s of sequences) {
    for (let i = 0; i + n <= s.length; i++) {
      const sub = s.substring(i, i + n);
      counts.set(sub, (counts.get(sub) ?? 0) + 1);
    }
  }
  return counts;
}

function countMatches(sequence: string, other: Map<string, number>, n: number): number {
  let total = 0;
  for (let i = 0; i + n <= sequence.length; i++)
    total += other.get(sequence.substring(i, i + n)) ?? 0;
  return total;
}

// Canvas colors encode the data (a common nucleotide scheme), so they are not theme tokens.
const NUCLEOTIDE_COLORS: {[n: string]: string} = {A: '#2ca02c', C: '#1f77b4', G: '#ff7f0e', T: '#d62728'};
const OTHER_COLOR = '#999999';
const CELL_PADDING = 2;
const LINE_HEIGHT = 13;

/** How many characters fit per line and how many lines to draw, filling lines in order. */
export function nucleotideLayout(length: number, w: number, h: number, charWidth: number):
  {perLine: number, lines: number} {
  const perLine = Math.max(1, Math.floor((w - 2 * CELL_PADDING) / charWidth));
  const maxLines = Math.max(1, Math.floor((h - 2 * CELL_PADDING) / LINE_HEIGHT));
  return {perLine, lines: Math.min(maxLines, Math.max(1, Math.ceil(length / perLine)))};
}

export class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name(): string { return 'Nucleotide cell renderer'; }

  get cellType(): string { return 'dna_nucleotide'; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    const text: string = gridCell.cell.value ?? '';
    g.font = '11px courier';
    g.textBaseline = 'middle';
    g.textAlign = 'left';
    const charWidth = g.measureText('A').width;
    const prefixLength = text.match(FASTA_PREFIX)?.[0].length ?? 0;
    const {perLine, lines} = nucleotideLayout(text.length, w, h, charWidth);
    // The block of lines is centered vertically, so a single line sits in the middle as before.
    const top = y + (h - lines * LINE_HEIGHT) / 2 + LINE_HEIGHT / 2;
    // Characters beyond the last line that fits are not drawn.
    const visible = Math.min(text.length, perLine * lines);
    for (let i = 0; i < visible; i++) {
      const n = text[i];
      g.fillStyle = i < prefixLength ? OTHER_COLOR : NUCLEOTIDE_COLORS[n.toUpperCase()] ?? OTHER_COLOR;
      g.fillText(n, x + CELL_PADDING + (i % perLine) * charWidth, top + Math.floor(i / perLine) * LINE_HEIGHT);
    }
  }
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_nucleotide
//meta.role: cellRenderer
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer(): DG.GridCellRenderer {
  return new NucleotideBoxCellRenderer();
}

//name: enaSequence
//tags: panel, widgets
//input: string cellText {semType: EnaID}
//output: widget result
//condition: true
export async function enaSequence(cellText: string): Promise<DG.Widget> {
  const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${encodeURIComponent(cellText)}`;
  const response = await grok.dapi.fetchProxy(url);
  if (!response.ok)
    return new DG.Widget(ui.divText(`ENA returned ${response.status} for ${cellText}`));
  const [header, ...lines] = (await response.text()).trim().split('\n');
  const sequence = ui.input.textArea('', {value: lines.join('')});
  return new DG.Widget(ui.box(ui.splitV([ui.divText(header), sequence.root])));
}

// The exercise suggests the embl text search, but its CON entries carry no SQ block, only a join of other
// entries; the fasta flavor of the same search returns the assembled sequence.
export async function _fetchENASequence(query: string, limit: number, length: number = 60): Promise<DG.DataFrame> {
  const params = new URLSearchParams({result: 'sequence', query, limit: `${limit}`});
  const response = await grok.dapi.fetchProxy(`https://www.ebi.ac.uk/ena/browser/api/fasta/textsearch?${params}`);
  if (!response.ok)
    throw new Error(`ENA search returned ${response.status}`);
  const ids: string[] = [];
  const sequences: string[] = [];
  for (const entry of (await response.text()).split('>').filter((e) => e.trim())) {
    const [header, ...lines] = entry.trim().split('\n');
    ids.push(header.split('|')[1]);
    sequences.push(lines.join('').substring(0, length));
  }
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ids),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', sequences),
  ]);
  df.col('ID')!.semType = 'EnaID';
  df.col('Sequence')!.semType = 'dna_nucleotide';
  return df;
}

//name: formENADataTable
//description: Searches ENA and creates a table of IDs and sequences
export function formENADataTable(): void {
  const queryInput = ui.input.string('Query', {value: 'coronavirus'});
  const lengthInput = ui.input.int('Sequence length', {value: 60});
  const rowsInput = ui.input.int('How many rows', {value: 100});
  const preview = ui.div();
  const previewButton = ui.button('Preview', async () => {
    ui.setUpdateIndicator(preview, true);
    try {
      const df = await _fetchENASequence(queryInput.value, 10, lengthInput.value!);
      preview.replaceChildren(DG.Viewer.grid(df).root);
    } finally {
      ui.setUpdateIndicator(preview, false);
    }
  });

  ui.dialog('Create sequences table')
    .add(ui.inputs([queryInput, lengthInput]))
    .add(previewButton)
    .add(preview)
    .add(ui.inputs([rowsInput]))
    .onOK(async () => {
      const df = await _fetchENASequence(queryInput.value, rowsInput.value!, lengthInput.value!);
      df.name = queryInput.value;
      grok.shell.addTableView(df);
    })
    .show();
}

//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: int N
//output: dataframe result
export function fuzzyJoin(df1: DG.DataFrame, df2: DG.DataFrame, N: number): DG.DataFrame {
  const col1 = df1.columns.bySemType('dna_nucleotide')!;
  const col2 = df2.columns.bySemType('dna_nucleotide')!;
  const df2Renamed = df2.clone();
  df2Renamed.col(col2.name)!.name = col1.name;
  const df = df1.append(df2Renamed);

  const seqs1 = col1.toList().map(normalize);
  const seqs2 = col2.toList().map(normalize);
  const subs1 = countSubsequences(seqs1, N);
  const subs2 = countSubsequences(seqs2, N);
  const counts = df.columns.addNew('Counts', DG.TYPE.INT);
  for (let i = 0; i < df.rowCount; i++) {
    counts.set(i, i < df1.rowCount ?
      countMatches(seqs1[i], subs2, N) :
      countMatches(seqs2[i - df1.rowCount], subs1, N));
  }

  grok.shell.addTableView(df);
  return df;
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(nucleotides: string): DG.Widget {
  return new DG.Widget(ui.divText(complement(nucleotides)));
}

// Named differently from the script, since function names are case-insensitive within a package namespace.
//name: countSubsequencePythonWrapper
//input: string sequence
//input: string subsequence
//output: int count
export async function countSubsequencePythonWrapper(sequence: string, subsequence: string): Promise<number> {
  return await grok.functions.call('GuillermoSequence:CountSubsequencePython', {sequence, subsequence});
}

//name: countSubsequenceTableAugment
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
export async function countSubsequenceTableAugment(
  sequences: DG.DataFrame, columnName: DG.Column, subsequence: string): Promise<void> {
  const df: DG.DataFrame = await grok.functions.call('GuillermoSequence:CountSubsequencePythonDataframe',
    {sequences, columnName, subsequence});
  const countCol = df.columns.byIndex(0);
  countCol.name = sequences.columns.getUnusedName(`N(${subsequence})`);
  sequences.columns.add(countCol);
}

//name: getOrders
//output: dataframe df
export async function getOrders(): Promise<DG.DataFrame> {
  return await grok.data.query('GuillermoSequence:ordersByCountry', {country: 'USA'});
}

//name: openTable1
//description: Opens a demo file with grok.data.getDemoTable
//input: string filepath = "bio/sars-cov-2.csv"
//output: dataframe df
export async function openTable1(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.getDemoTable(filepath);
  grok.shell.addTableView(df);
  return df;
}

//name: openTable2
//description: Opens a demo file with grok.data.files.openTable
//input: string filepath = "bio/sars-cov-2.csv"
//output: dataframe df
export async function openTable2(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.files.openTable(`System:DemoFiles/${filepath}`);
  grok.shell.addTableView(df);
  return df;
}

//name: openTable3
//description: Opens a demo file with the OpenServerFile command
//input: string filepath = "bio/sars-cov-2.csv"
//output: dataframe df
export async function openTable3(filepath: string): Promise<DG.DataFrame> {
  const df: DG.DataFrame = (await grok.functions.eval(`OpenServerFile("System:DemoFiles/${filepath}")`))[0];
  grok.shell.addTableView(df);
  return df;
}

//name: addTables
//description: Opens every CSV file shipped in the package files folder
export async function addTables(): Promise<void> {
  const files = await _package.files.list('', true);
  for (const file of files.filter((f) => f.extension === 'csv'))
    grok.shell.addTableView(await _package.files.readCsv(file));
}
