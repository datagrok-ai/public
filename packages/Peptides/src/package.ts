/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {splitAlignedPeptides} from './splitAligned';
import {describe} from './SARViewer/describe';

export const _package = new DG.Package();


// function main(tableName, localTables) {
//   grok.shell.info(tableName);
//   let peptides = localTables.find(({name}) => name === tableName);
//   let view = grok.shell.addTableView(peptides);
//   grok.shell.v = view;
// }

async function main() {
  const pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides = await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562')
  const path = _package.webRoot + 'files/' + 'aligned.csv';
  const peptides = (await grok.data.loadTable(path));


  // let names = [
  //   "N", "a3", "a4", "a5", "a6",
  //   "a7", "a8", "a9", "a10",
  //   "a11", "a12", "a13", "a14",
  //   "a15", "a16", "a17", "C", "Activity"
  // ];
  // let htmlNames = [
  //   "a3", "a4", "a5", "a6",
  //   "a7", "a8", "a9", "a10",
  //   "a11", "a12", "a13", "a14",
  //   "a15", "a16", "a17"
  // ];

  const view = grok.shell.addTableView(peptides);

  const grid = await describe(peptides, 'Activity');

  if (grid !== null) {
    view.addViewer(grid);
  }

  // for(let i = 0; i < names.length; i++){
  //   if(htmlNames.includes(names[i])){
  //     view.grid.columns.byName(names[i]).width = 50;
  //     view.grid.columns.byName(names[i]).cellType = 'html';
  //   }
  // }

  // let css = {
  //   style:{
  //   textAlign:'center',
  //   fontWeight:'bold',
  //   marginTop:'5px'
  //   }
  // };

  // let css2 = {
  //   style:{
  //   textAlign:'center',
  //   marginTop:'5px',
  //   color: "red"
  //   }
  // };

  // let css3 = {
  //   style:{
  //   textAlign:'center',
  //   marginTop:'5px',
  //   color: "green"
  //   }
  // };

  // let aaRgx = new RegExp("^[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}$");
  // let adRgx = new RegExp("^d[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}$");
  // let amRgx = new RegExp("^[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}\\(.+\\)$");
  // let admRgx = new RegExp("^d[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}\\(.+\\)$$");
  // let aNameRgx = new RegExp("^a\\d.*$");

  // view.grid.onCellPrepare(function (gc) {
  //   if (gc.isTableCell && gc.gridColumn.name === 'N'){
  //     gc.style.backColor = 0xff1f77b4;
  //   }
  //   if (gc.isTableCell && gc.gridColumn.name === 'C'){
  //     gc.style.backColor = 0xffff7f00;
  //   }
  //   if(gc.isTableCell && !aaRgx.test(gc.cell.value) && !amRgx.test(gc.cell.value) && aNameRgx.test(gc.gridColumn.name)){
  //     gc.style.element = ui.divText(gc.cell.value, css2);
  //   }
  //   if(gc.isTableCell && aaRgx.test(gc.cell.value)){
  //     console.log(gc)
  //     gc.style.element = ui.divText(gc.cell.value, css);
  //   }
  //   if(gc.isTableCell && amRgx.test(gc.cell.value)){
  //     console.log(gc)
  //     gc.style.element = ui.divText(gc.cell.value, css);
  //   }
  //   if(gc.isTableCell && adRgx.test(gc.cell.value)){
  //     console.log(gc)
  //     gc.style.element = ui.divText(gc.cell.value, css3);
  //   }
  //   if(gc.isTableCell && admRgx.test(gc.cell.value)){
  //     console.log(gc)
  //     gc.style.element = ui.divText(gc.cell.value, css3);
  //   }

  // });

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {
  const appDescription = ui.info(
    [
      ui.span(['For more deatails see LINK ']),
      ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
      ui.divText('Select the corresponding .csv table with peptide sequences'),
    ], 'Transform peptide sequence data to research insights'
  );
  const annotationViewerDiv = ui.div();

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  // const localTables = grok.shell.tables;
  // const namesOfLocalTables = localTables.map((df) => df.name);

  // let chosenFile = ui.choiceInput('File', '', namesOfLocalTables, () => {
  //   main(chosenFile.stringValue, localTables);
  // });

  const chosenFile = ui.button('Open demo', ()=> main());

  const mainDiv = ui.div();
  grok.shell.newView('Peptides', [
    appDescription,
    ui.h2('Choose .csv file'),
    ui.div([
      ui.block25([
        ui.inputs([
          chosenFile,
        ]),
      ]),
      ui.block75([annotationViewerDiv]),
    ]),
    mainDiv,
  ]);

  //let table = DG.DataFrame.create();
  //table.name = 'Peptides';
  //let view = grok.shell.addTableView(table);
}


//name: peptidesSAR
//tags: demo, peptides
//description:
//top-menu: Chem | peptides
//input: dataframe df {semType: alignedSequence}
export function peptidesSAR(df: DG.DataFrame) {
  grok.shell.warning(df.getCol('Activity').toList().length.toString());

  //let colSmiles = DG.Column.fromList("string", "smiles", smiles.toList());
  //let dfSmiles = DG.DataFrame.fromColumns([colSmiles]);


  //df.columns.insert(x_coord);
  // df.columns.insert(y_coord);
  // df.columns.insert(sali);

  // df.columns.byName('x_coord').name = '~x_coord';
  // df.columns.byName('y_coord').name = '~y_coord';

  // let view = grok.shell.getTableView(df.name);
  // let sp = view.addViewer(DG.Viewer.scatterPlot(df, {
  //   xColumnName: '~x_coord',
  //   yColumnName: '~y_coord',
  //   size: 'sali'
  // }));

  // sp.props.showXSelector = false;
  // sp.props.showYSelector = false;
  // sp.props.showSizeSelector = false;
  // sp.props.showColorSelector = false;
  // sp.props.markerMinSize = 5;
  // sp.props.markerMaxSize = 25;
  // sp.props.colorColumnName = 'activity';
}

//name: Split Sequence
//input: column peptideColumn {semType: alignedSequence}
//tags: panel
//output: widget result
export function splitAlignedSequence(peptideColumn: DG.Column) {
  grok.shell.addTableView(splitAlignedPeptides(peptideColumn));
}
