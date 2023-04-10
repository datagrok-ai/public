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
//input: string nucleotides {semType: dna_sequences}
//output: string result {semType: dna_sequences}
export function complement(nucleotides: string): string {
  let out = '';

  const getComplement = (nuc:string) =>{
    switch (nuc.toLowerCase()) {
    case 'a':
      return 'T';
    case 't':
      return 'A';
    case 'g':
      return 'C';
    case 'c':
      return 'G';
    case ' ':
      return ' ';
    default:
      throw Error('DNA nucleotydes must be either A/C/G/T');
    }
  };
  for (let i = 0; i<nucleotides.length; i++)
    out += getComplement(nucleotides[i]);

  return out;
}

//name: complement Sequence
//tags: panel, widgets
//input: string nucleotides {semType: dna_sequences}
//output: widget result
//condition: true
export function complementWidget(nucleotides:string) {
  return nucleotides ? new DG.Widget(ui.divText(complement(nucleotides))) : new DG.Widget(ui.divText('value is empty'));
}


//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: int N
export function fuzzyJoin(df1: DG.DataFrame, df2: DG.DataFrame, N: number) {
  const col1 = df1.columns.bySemTypesExact(['dna_sequences'])?.[0];
  const col2 = df2.columns.bySemTypesExact(['dna_sequences'])?.[0];
  if (!col1 || !col2)
    throw new EvalError('one of the dataframes does not contain dna_sequences semantic type');

  const df = df1.append(df2);
  const newCol = df.columns.addNew('Counts', 'int');
  for (let i = 0; i< (col1?.length?? 0); i++) {
    let c = 0;

    for (let j = 0; j < col1.get(i).length - N; j++) {
      const subsequence = col1.get(i).substring(j, j+N);
      let count = 0;
      for (let k = 0; k < col2.length; k++)
        count +=countAllOccurances(col2.get(k), subsequence);

      c+=count;
    }
    newCol.set(i, c);
  }

  grok.shell.addTableView(df);
}

function countAllOccurances(sequence:any, subsequence:any) {
  let count = 0;
  for (let i=0; i<=sequence.length - subsequence.length; i++) {
    if (sequence.startsWith(subsequence, i))
      count +=1;
  }
  return count;
}


export class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name() {return 'Nucleotide cell renderer';}
  get cellType() {return 'dna_nucleotide';}
  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    const cellProps = {
      cellWidth: 300,
      fontSize: 11,
    };


    const seq = gridCell.cell.value;
    const ctx = g.canvas.getContext('2d')!;
    ctx.font = `${cellProps.fontSize}px courier bold`;
    ctx.fillStyle = '#FF0000';
    ctx.textAlign = 'center';

    //const rowsToRender = Math.ceil((seq.length??0) * cellProps.fontSize / cellProps.cellWidth);
    //const rowHeight = cellProps.fontSize * 1.5 * rowsToRender;
    const maxWidth = Math.max(gridCell.gridColumn.width, cellProps.cellWidth) + 20;
    //const maxHeight = rowHeight;

    gridCell.gridColumn.width = maxWidth;
    //ctx.canvas.height = maxHeight;
    for (let i = 0; i < gridCell.cell.value.length; i++) {
      const nuc = seq[i].toUpperCase();
      switch (nuc) {
      case 'A':
        ctx.fillStyle = '#00FF00';
        break;
      case 'C':
        ctx.fillStyle = '#0000FF';
        break;
      case 'T':
        ctx.fillStyle = '#FF0000';
        break;
      default:
        ctx.fillStyle = '#000000';
      }
      ctx.fillText(nuc, x+(i%30+1)*cellProps.fontSize*0.9,
        y+cellProps.fontSize*1.5 + cellProps.fontSize*1.5 * Math.floor(i/30),
        cellProps.cellWidth);
    }
  }
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_sequences
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer() {
  return new NucleotideBoxCellRenderer();
}


//name: ENA Sequence
//tags: panel, widgets
//input: string cellText {semType: EnaID}
//output: widget result
//condition: true
export async function enaSequence(cellText: string) {
  const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${cellText}`;
  const fasta = await (await grok.dapi.fetchProxy(url)).text();
  if (!fasta || !fasta.length || !fasta.startsWith('>'))
    throw Error('incorrect fasta format');

  const data = fasta.split('\n');
  const name = data[0];
  const sequence = data.reduce((acc, cur, j) => {
    if (j!==0)
      return acc + cur;
    return acc;
  }, '');

  //const [name, sequence] = fasta.split('\n');


  return new DG.Widget(ui.box(
    ui.splitV([ui.divText(name.substring(1, name.length)), ui.textInput('', sequence).root]),

  ));
}

interface SequenceInfo{
  name:string;
  sequence:string;
}

export async function getByEnas(enas:string[]) {
  const out:SequenceInfo[] = [];
  for (let i = 0; i<enas.length; i++) {
    const ena = enas[i];
    const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${ena}`;
    const fasta = await(await grok.dapi.fetchProxy(url)).text();
    if (!fasta || !fasta.length || !fasta.startsWith('>'))
      throw Error('incorrect fasta format');

    const sequence = fasta.split('\n').reduce((cur, acc, j) => {
      if (j!==0)
        return acc + cur;
      return acc;
    }, '');
    out.push({name: ena, sequence});
  }
  return out;
}

//name: Form ENA Table by the given topic
export async function formENADataTable() {
  const getUrl = (topic:string = 'coronavirus', limit:number = 10) =>
    `https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=${topic}&limit=${limit}`;

  const onPreview = async (grid:DG.Grid, topic?:string, limit?:number) =>{
    const enas = await getByTopic(topic, limit);
    console.log(enas);
    const sequences = await getByEnas(enas);
    console.log(sequences);
    const sequenceCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', sequences.map((seq) => seq.sequence));
    sequenceCol.semType = 'dna_sequences';
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', sequences.map((seq) => seq.name)),
      sequenceCol,
    ]);
    grid.dataFrame = df;
    return df;
  };

  const getByTopic = async (topic?:string, limit?:number) => {
    const url = getUrl(topic, limit);
    const data = await (await grok.dapi.fetchProxy(url)).text();
    const enas:string[] = [];
    data.split('\n').forEach((line) =>{
      if (line.startsWith('ID'))
        enas.push(line.substring(5, 13));
    });
    return enas;
  };

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', []),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', []),
  ]);

  const grid = DG.Viewer.grid(df);
  const limitInput = ui.intInput('How many rows: ', 100);
  const queryInput = ui.stringInput('Query: ', 'coronavirus');
  const button = ui.button('Preview', () => onPreview(grid));
  ui.dialog('Create sequences table')
    .add(ui.splitV([

      queryInput.root,
      ui.splitH([
        ui.div([limitInput]),
        button,
      ]),


      ui.div([grid]),

    ]))
    .onOK(async () => {
      const n = limitInput.value ?? 100;
      const topic = queryInput.value ?? 'coronavirus';
      const df = await onPreview(grid, topic, n);
      grok.shell.addTableView(df);
    })
    .show();
}

