/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as nu from './nucleotide-utils';
import * as ena from './ena';
import {Nucleotide} from './types';

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

  // const subsequenceLists = subsequencesCol.toList().map((seq: string): string[] => (
  //   nu.generateSubsequences(seq.replaceAll(/\s+/g, ''), N)
  // ));

  // for (let i = 0; i < subsequenceLists.length; ++i) {
  //   const isFirstfHalf = i < df1Size;
  //   const jEnd = isFirstfHalf ? subsequenceLists.length : df1Size;
  //   let matches = 0;
  //   for (let j = isFirstfHalf ? df1Size : 0; j < jEnd; ++j) {
  //     for (const seqEl of subsequenceLists[i])
  //       matches += subsequenceLists[j].filter((n: string): boolean => n === seqEl).length;
  //   }
  //   countCol.set(i, matches);
  // }

  const subsequenceLists = subsequencesCol.toList().map((seq: string): Map<string, number> => (
    nu.countSubsequences(seq.replaceAll(/\s+/g, ''), N)
  ));

  for (let i = 0; i < subsequenceLists.length; ++i) {
    const isFirstfHalf = i < df1Size;
    const jEnd = isFirstfHalf ? subsequenceLists.length : df1Size;
    let matches = 0;
    for (let j = isFirstfHalf ? df1Size : 0; j < jEnd; ++j) {
      for (const [seqEl, seqCount] of subsequenceLists[i])
        matches += seqCount * (subsequenceLists[j].get(seqEl) ?? 0);
    }
    countCol.set(i, matches);
  }

  // grok.shell.addTableView(df);
  return df;
}

export class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name() {return 'Nucleotide cell renderer';}
  get cellType() {return nu.NUCLEOTIDE_SEMTYPE;}
  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    const sequence: string[] = gridCell.cell.value;
    const ctx = g.canvas.getContext('2d');
    if (!ctx) return;

    // ctx.font = '11px courier';
    ctx.font = cellStyle.font;

    const paddingX = 3;
    const paddingY = 3;

    const charOffsetX = 8;
    const alternativeCharOffsetX = 6;
    const charOffsetY = 13;

    const x2 = x + w - paddingX;
    const y2 = y - paddingY;

    let curX = x + paddingX;
    let curY = y - h * .7 + paddingY;
    for (const nucleotideRaw of sequence) {
      if (curY > y2) break;
      const nucleotide = nucleotideRaw.toUpperCase() as Nucleotide;

      ctx.fillStyle = nu.NUCLEOTIDE_COLORS[nucleotide] ?? 'black';
      ctx.fillText(nucleotideRaw, curX, curY);

      const isAtrernativeOffset = nucleotide === 'T';
      curX += isAtrernativeOffset ? alternativeCharOffsetX : charOffsetX;

      if (curX >= x2) {
        curX = x + paddingX;
        curY += charOffsetY;
      }
    }
  }
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
  const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${cellText}`;
  const fasta = await grok.dapi.fetchProxy(url)
    .then((res: Response): Promise<string> => res.text());
  const enaSequence = ena.parseENASequenceFasta(fasta);
  // console.log('enaSequence:', enaSequence);

  return new DG.Widget(ui.box(
    // ui.splitV([
    ui.divV([
      ui.h1(enaSequence?.id ?? 'No ENA id'),
      ui.divText(enaSequence?.sequence ?? 'Sequence is empty'),

      // ui.divV([ui.h3('code'), ui.divText(enaSequence?.code ?? 'is empty')]),
      // ui.divV([ui.h3('extra'), ui.divText(enaSequence?.extra ?? 'is empty')]),
      // ui.divV([ui.h3('description'), ui.divText(enaSequence?.description ?? 'is empty')]),
      // ui.divV([ui.h3('genBank'), ui.divText(enaSequence?.genBank ?? 'is empty')]),
      // ui.divV([ui.h3('sequence'), ui.input.textArea(enaSequence?.sequence ?? 'is empty')]),

      // ui.divV([ui.h3('sequence'), ui.divText(enaSequence?.sequence ?? 'is empty')]),
    ], {style: {gap: '0'}}),
  ));
}

//name: _fetchENASequence
//input: string query
//input: int limit = 3
//input: int offset = 0
//output: dataframe result
export async function _fetchENASequence(query: string, limit: number, offset: number): Promise<DG.DataFrame> {
  const queryParams = `result=sequence&query=${query}&limit=${limit}&offset=${offset}`;
  const idsUrl = `https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?${queryParams}`;
  const seqsUrl =`https://www.ebi.ac.uk/ena/browser/api/fasta/textsearch?${queryParams}`;

  const idsRequest: Promise<string[]> = grok.dapi.fetchProxy(idsUrl)
    .then((res: Response): Promise<string> => res.text())
    .then(ena.searchENAIdsEmbs);

  const seqsReqest: Promise<string[]> = grok.dapi.fetchProxy(seqsUrl)
    .then((res: Response): Promise<string> => res.text())
    .then(ena.parseENAMultipleSequencesFasta);

  const [ids, seqs] = await Promise.all([idsRequest, seqsReqest]);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ids),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', seqs),
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

  const grid = DG.Viewer.grid(previewDf);
  const limitInput = ui.input.int('How many rows: ', {value: 100});
  const queryInput = ui.input.string('Query: ', {value: 'coronavirus'});
  const offsetInput = ui.input.int('Sequence offset: ', {value: 0});

  const createDf = (): Promise<DG.DataFrame> => (
    _fetchENASequence(queryInput.value, limitInput.value, offsetInput.value)
  );

  const button = ui.button('Preview', async (): Promise<void> => {
    if (previewDf.rowCount > 0)
      previewDf.rows.removeAt(0, previewDf.rowCount);
    previewDf.append(await createDf(), true);
  });

  ui.dialog('Create sequences table')
    .add(ui.splitV([
      ui.splitH([
        ui.span([queryInput.root]),
        button,
      ]),
      ui.div([grid]),
      ui.div([limitInput, offsetInput]),
    ]))
    .onOK(async (): Promise<void> => {
      const df = previewDf.rowCount > 0 ? previewDf : await createDf();
      grok.shell.addTableView(df);
    })
    .show();
}
