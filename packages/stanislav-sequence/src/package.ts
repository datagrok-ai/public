/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { swapTwoWays } from './swap';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: complement
//input: string nucleotides
//output: string result
export async function complement(nucleotide: string): Promise<string> {
  let result = swapTwoWays(nucleotide, 'A', 'T');
  result = swapTwoWays(nucleotide, 'C', 'G');
  return result;
}

//name: dna_nucleotide 
//input: string nucleotides {semType: dna_nucleotide}
//output: string result {semType: dna_nucleotide}
export function dna_nucleotide(nucleotide: string): string {
  // your code goes here
  return nucleotide;
}

//name: complementWidget
//tags: panel, widgets
//input: string nucleotides {semType: dna_nucleotide}
//output: widget result
//condition: true
export function complementWidget(nucleotide: any): any {
  // your code goes here
  return new DG.Widget(ui.divText(dna_nucleotide(nucleotide)));;
}

//name: executeFunction
//tags: panel, widgets
//input: string sequence 
//input: string subsequence    
//output: int result 
export async function executeCountSubsequence(sequence: string, subsequence: string): Promise<number> {
  return await grok.functions.call('StanislavSequence:CountSubsequencePython', { sequence, subsequence });
}

//name: getOrders
//input: string country 
//output: dataframe df
export async function getOrders(country: string): Promise<DG.DataFrame> {
  return await grok.data.query(`StanislavSequence:OrdersByCountry`, { country: country });
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
  const df = await grok.data.files.openTable(`System:/${filepath}`);
  grok.shell.addTableView(df);
  return df;
}

//input: string filepath
//output: dataframe df
export async function openTable3(filepath: string): Promise<DG.DataFrame> {
  const df = (await (grok.functions.eval(`OpenServerFile("System:DemoFiles/${filepath}")`)))[0];
  grok.shell.addTableView(df);
  return df;
}

//input: string filepath
//output: dataframe df
export async function openTable4(filepath: string): Promise<DG.DataFrame> {
  const df = await grok.data.loadTable(`https://dev.datagrok.ai/${filepath}`);
  grok.shell.addTableView(df);
  return df;
}

//name: Add Tables
//output: string t
export async function addTables(): Promise<String> {
  const files = await _package.files.list('', true);

  const csvFiles = files.filter((f) => f.extension === 'csv');
  let t = "";
  // Load every table and add a view for it
  for (const file of csvFiles) {
    const df = await _package.files.readCsv(file.fileName);
    t = t + file + ' ';
    grok.shell.addTableView(df);
    // const df = await grok.data.loadTable(`${_package.webRoot}${file.path}`);
    // const df = await grok.data.files.openTable(`System:AppData/${_package.name}/${file.fileName}`);
  }
  return t;
}

//name: fuzzyJoin
//input: dataframe df1
//input: dataframe df2
//input: column col
//input: int N 
//output: dataframe result
export function fuzzyJoin(df1: DG.DataFrame, df2: DG.DataFrame, col: DG.Column, N: number) {
  const col1 = df1.columns.byName(col.name).toList();
  const col2 = df2.columns.byName(col.name).toList();
  const col1Length = col1.length;

  let resultDF = df1.append(df2);
  resultDF.columns.addNew("Subsequence Count", "int")

  let i = 0;
  for (let item of resultDF.rows) {
    if (i < col1Length) {
      item["Subsequence Count"] = getCountOfSubsequencesInColumn(item[col.name], N, col2)
    }
    else {
      item["Subsequence Count"] = getCountOfSubsequencesInColumn(item[col.name], N, col1)
    }
    i++;
  }

  grok.shell.addTableView(resultDF);
  return resultDF;
}

function getCountOfSubsequencesInColumn(subsequenceRowString: string, subseqSize: number, colValues: any): number {
  let result = 0;
  let subseqItems: Set<string> = new Set<string>();

  for (let i = 0; i + subseqSize <= subsequenceRowString.length; i++) {
    subseqItems.add(subsequenceRowString.substring(i, i + subseqSize));
  }

  for (let item of colValues) {
    for (let subseqItem of subseqItems) {
      result += (item.match(subseqItem) || []).length;
    }
  }
  return result;
}

//name: nucleotideBoxCellRenderer
//tags: cellRenderer
//meta.cellType: dna_nucleotide
//output: grid_cell_renderer result
export function nucleotideBoxCellRenderer() {
  return new NucleotideBoxCellRenderer();
}

export class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Nucleotide cell renderer'; }
  get cellType() { return 'dna_nucleotide'; }


  getDefaultSize(gridColumn: DG.GridColumn): { width?: number | null, height?: number | null } {
    return { width: 200, height: 30 };
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    let seq = gridCell.cell.value;
    let ctx = g.canvas.getContext('2d')!;
    ctx.font = '11px courier';
    // ...
    let localx = x + 22
    for (let i = 0; i < gridCell.cell.value.length; i++) {
      localx += 11;
      switch (gridCell.cell.value[i]) {
        case 't':
        case 'T':
          ctx.fillStyle = "red";
          break;
        case 'a':
        case 'A':
          ctx.fillStyle = "blue";
          break;
        case 'c':
        case 'C':
          ctx.fillStyle = "green";
          break;

        case 'g':
        case 'G':
          ctx.fillStyle = "black";
          break;
      }
      ctx.fillText(gridCell.cell.value[i], localx, y, h);
    }

  }
}

//name: ENA Sequence
//tags: panel, widgets
//input: string cellText {semType: EnaID}
//output: widget result
//condition: true
export async function enaSequence(cellText: string) {
  const url = `https://www.ebi.ac.uk/ena/browser/api/fasta/${cellText}`;
  const fasta = await (await grok.dapi.fetchProxy(url)).text();
  return new DG.Widget(ui.box(
    ui.textInput(cellText, fasta),
    // ... the widget controls are composed here
  ));
}

//name: form ENA Data Table
//tags: panel, widgets 
//output: widget result
//condition: true
export async function formENADataTable() {

  let df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ["a", "b"]),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sequence', ["c", "d"])
  ]);

  let inputData = 'coronavirus'
  let grid = DG.Viewer.grid(df);

  let limitInput = ui.intInput('How many rows: ', 100);
  let queryInput = ui.stringInput('Query: ', inputData);
 
  let button = ui.button('Preview', () => { 
  
    browseDataForENAID(queryInput.value, limitInput.value|| 10).then((result) =>{
      df = DG.DataFrame.fromColumns( [ DG.Column.fromStrings("EnaId", (result)) ] );   
      grid.dataFrame = df
    });
  });
  
  ui.dialog('Create sequences table')
    .add(ui.splitV([
      ui.splitH([
        ui.span([queryInput.root]),
        button
      ]),
      grid.root,
      ui.div([limitInput])
    ]))
    .onOK(() => {
      browseDataForENAID(queryInput.value, limitInput.value|| 10).then((result) =>{
        df = DG.DataFrame.fromColumns( [ DG.Column.fromStrings("EnaId", (result)) ] );  
        grok.shell.addTableView(df);
      });
    })
    .show();
}

function fetchEnaIdSequence(queryString: string, elementsCount: number): string {
  return ""
}

async function browseDataForENAID(queryString: string, elementsCount: number): Promise<Array<string>> {
  let browseDataForENAIDResult = []
  const queryUrl = `https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=${queryString}&limit=${elementsCount}&annotationOnly=false` 
  const queryResult = await (await grok.dapi.fetchProxy(queryUrl)).text()
  const EnaIDRegex = "AC   ([A-Z]{2}[0-9]{6}(.[0-9])?)";
  console.log(queryResult);

  for(let queryResultLine of queryResult.split(/[\r\n]+/)){

    const matchResult = queryResultLine.match(EnaIDRegex) || [];
    if(matchResult.length  > 1 ){
      browseDataForENAIDResult.push(matchResult[1])
    }
  }
  return browseDataForENAIDResult
}


// export class TagsCellRenderer extends DG.GridCellRenderer {
//   get name() { return 'Tags'; }

//   get cellType() { return 'Tags'; }

//   render(
//     g: CanvasRenderingContext2D,
//     x: number, y: number, w: number, h: number,
//     gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
//   ) {
//     // somehow dart comes null here, will need to investigate it and fix it, now just a workaround
//     if (!gridCell.gridColumn.dart)
//       return;
//     const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());

//     let cx = 2;
//     let cy = 3;
//     g.textBaseline = 'top';

//     const getColor = (tag: string) => {
//       const colors = gridCell.gridColumn.temp['catColors'] ??= {};
//       if (colors[tag] || (colors[tag] === 0))
//         return colors[tag];

//       const keys = Object.keys(colors);
//       colors[tag] ??= DG.Color.getCategoricalColor(keys.length);
//     };

//     for (const tag of values) {
//       const width = g.measureText(tag).width;
//       const drawTag = () => {
//         const color = getColor(tag);
//         g.fillStyle = DG.Color.toHtml(color);
//         g.roundRect(x + cx, y + cy, width + 4, 16, 4);
//         g.fill();

//         g.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
//         g.fillText(tag, x + cx + 2, y + cy + 1);
//       };

//       if (cx + width <= w) {
//         drawTag();
//         cx += width + 8;
//       }
//       else {
//         cx = 2;
//         cy += 16;
//         drawTag();
//       }
//     }
//   }
// }
