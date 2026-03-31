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
    {sequence, subsequence}
  );
}

//name: CountSubsequenceTableAugment
//input: dataframe sequences
//input: column columnName
//input: string subsequence = "acc"
export async function countSubsequenceTableAugment(sequences: DG.DataFrame, columnName: DG.Column, subsequence: string): Promise<void> {
  const df = await grok.functions.call('valerij_developer_exercise_1:CountSubsequencePythonDataframe', {sequences, columnName, subsequence});
  const countCol = df.columns.byIndex(0);
  countCol.name = `N(${subsequence})`;
  sequences.columns.insert(countCol);
}

//name: getOrders
//output: dataframe df
export async function getOrders() {
  return await grok.data.query('valerij_developer_exercise_1:ordersByCountry', {country: 'USA'});
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

//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: int N
//output: dataframe result
export function fuzzyJoin(df1: DG.DataFrame, df2: DG.DataFrame, N: number): DG.DataFrame {
  const col1 = findFirstDnaColumn(df1);
  const col2 = findFirstDnaColumn(df2);

  if (col1 == null)
    throw new Error('df1 does not contain a column with semType dna_nucleotide');

  if (col2 == null)
    throw new Error('df2 does not contain a column with semType dna_nucleotide');

  if (N <= 0)
    throw new Error('N must be greater than 0');

  const seqColName = col1.name;

  const left = df1.clone();
  const right = df2.clone();

  const rightSeqCol = right.col(col2.name);
  if (rightSeqCol == null)
    throw new Error(`Column ${col2.name} was not found in df2 clone`);

  rightSeqCol.name = seqColName;
  rightSeqCol.semType = 'dna_nucleotide';

  const result = left.clone();
  result.append(right, true);

  const leftSequences = getNormalizedSequences(df1, col1.name);
  const rightSequences = getNormalizedSequences(df2, col2.name);

  const countsCol = result.columns.addNewInt('Counts');

  for (let i = 0; i < df1.rowCount; i++) {
    const source = leftSequences[i];
    countsCol.set(i, fuzzyCount(source, rightSequences, N));
  }

  for (let i = 0; i < df2.rowCount; i++) {
    const source = rightSequences[i];
    countsCol.set(df1.rowCount + i, fuzzyCount(source, leftSequences, N));
  }

  grok.shell.addTableView(result);
  return result;
}

function findFirstDnaColumn(df: DG.DataFrame): DG.Column | null {
  for (let i = 0; i < df.columns.length; i++) {
    const col = df.columns.byIndex(i);
    if (col.semType === 'dna_nucleotide')
      return col;
  }

  const sequenceCol = df.col('sequence');
  if (sequenceCol != null) {
    sequenceCol.semType = 'dna_nucleotide';
    return sequenceCol;
  }

  return null;
}

function getNormalizedSequences(df: DG.DataFrame, colName: string): string[] {
  const col = df.col(colName);
  if (col == null)
    throw new Error(`Column ${colName} not found`);

  const values: string[] = [];
  for (let i = 0; i < df.rowCount; i++)
    values.push(normalizeSequence(String(col.get(i) ?? '')));

  return values;
}


function normalizeSequence(value: string): string {
  return value
    .toUpperCase()
    .replace(/^FASTA:\s*/, '')
    .trim();
}

function getSubsequences(sequence: string, n: number): string[] {
  if (sequence.length < n)
    return [];

  const unique = new Set<string>();
  for (let i = 0; i <= sequence.length - n; i++)
    unique.add(sequence.substring(i, i + n));

  return Array.from(unique);
}

function countOccurrences(sequence: string, subsequence: string): number {
  let count = 0;
  for (let i = 0; i <= sequence.length - subsequence.length; i++) {
    if (sequence.substring(i, i + subsequence.length) === subsequence)
      count++;
  }
  return count;
}

function fuzzyCount(source: string, targets: string[], n: number): number {
  const subsequences = getSubsequences(source, n);
  let total = 0;

  for (const subsequence of subsequences) {
    for (const target of targets)
      total += countOccurrences(target, subsequence);
  }

  return total;
}

export class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Nucleotide cell renderer'; }
  get cellType() { return 'dna_nucleotide'; }

  render(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle
  ) {
    const raw = String(gridCell.cell.value ?? '')
      .replace(/^FASTA:\s*/i, '')
      .trim();

    const ctx = g.canvas.getContext('2d');
    if (!ctx)
      return;

    ctx.save();

    ctx.beginPath();
    ctx.rect(x, y, w, h);
    ctx.clip();

    ctx.fillStyle = 'white';
    ctx.fillRect(x, y, w, h);

    const fontSize = Math.max(10, Math.min(14, h - 4));
    ctx.font = `${fontSize}px monospace`;
    ctx.textBaseline = 'top';

    const lineHeight = fontSize + 2;
    const startX = x + 4;
    let cursorX = startX;
    let cursorY = y + 4;

    for (let i = 0; i < raw.length; i++) {
      const ch = raw[i];

      if (ch === '\n') {
        cursorX = startX;
        cursorY += lineHeight;
        continue;
      }

      const width = ctx.measureText(ch).width;

      if (cursorX + width > x + w - 4) {
        cursorX = startX;
        cursorY += lineHeight;
      }

      ctx.fillStyle = nucleotideColor(ch);
      ctx.fillText(ch, cursorX, cursorY);
      cursorX += width;
    }

    ctx.restore();
  }
}

function nucleotideColor(ch: string): string {
  switch (ch.toUpperCase()) {
  case 'A': return '#4CAF50'; // green
  case 'T': return '#E53935'; // red
  case 'C': return '#1E88E5'; // blue
  case 'G': return '#222222'; // black/dark gray
  default: return '#666666';
  }
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_nucleotide
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer() {
  return new NucleotideBoxCellRenderer();
}
