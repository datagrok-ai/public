/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

// function main(tableName, localTables) {
//   grok.shell.info(tableName);
//   let peptides = localTables.find(({name}) => name === tableName);
//   let view = grok.shell.addTableView(peptides);
//   grok.shell.v = view;
// }

async function main() {
  let pi = DG.TaskBarProgressIndicator.create('Loading Peptides');
  //let peptides = await grok.data.loadTable('https://datagrok.jnj.com/p/ejaeger.il23peptideidp5562/il-23_peptide_idp-5562')
  let path = _package.webRoot + 'files/' + 'data.csv';

  let peptidesRaw = (await grok.data.loadTable(path))
  
  let names = [
    "N", "a3", "a4", "a5", "a6",
    "a7", "a8", "a9", "a10",
    "a11", "a12", "a13", "a14",
    "a15", "a16", "a17", "C", "Activity"
  ];
  let htmlNames = [
    "a3", "a4", "a5", "a6",
    "a7", "a8", "a9", "a10",
    "a11", "a12", "a13", "a14",
    "a15", "a16", "a17"
  ];

  let peptides = DG.DataFrame.fromColumns(peptidesRaw.columns.byNames(names));
  
  let view = grok.shell.addTableView(peptides);

  for(let i = 0; i < names.length; i++){
    if(htmlNames.includes(names[i])){
      view.grid.columns.byName(names[i]).width = 50;
      view.grid.columns.byName(names[i]).cellType = 'html';
    }
  }
  // view.grid.columns.byName('N').width = 50;
  // view.grid.columns.byName('a3').width = 50;
  // view.grid.columns.byName('a3').cellType = 'html';
  // view.grid.columns.byName('a4').width = 50;
  // view.grid.columns.byName('a4').cellType = 'html';
  // view.grid.columns.byName('a5').width = 50;
  // view.grid.columns.byName('a5').cellType = 'html';
  // view.grid.columns.byName('a6').width = 50;
  // view.grid.columns.byName('a6').cellType = 'html';
  // view.grid.columns.byName('a7').width = 50;
  // view.grid.columns.byName('a7').cellType = 'html';
  // view.grid.columns.byName('a8').width = 50;
  // view.grid.columns.byName('a8').cellType = 'html';
  // view.grid.columns.byName('a9').width = 50;
  // view.grid.columns.byName('a9').cellType = 'html';
  // view.grid.columns.byName('a10').width = 50;
  // view.grid.columns.byName('a10').cellType = 'html';
  // view.grid.columns.byName('a11').width = 50;
  // view.grid.columns.byName('a11').cellType = 'html';
  // view.grid.columns.byName('a12').width = 50;
  // view.grid.columns.byName('a12').cellType = 'html';
  // view.grid.columns.byName('a13').width = 50;
  // view.grid.columns.byName('a13').cellType = 'html';
  // view.grid.columns.byName('a14').width = 50;
  // view.grid.columns.byName('a14').cellType = 'html';
  // view.grid.columns.byName('a15').width = 50;
  // view.grid.columns.byName('a15').cellType = 'html';
  // view.grid.columns.byName('a16').width = 50;
  // view.grid.columns.byName('a16').cellType = 'html';
  // view.grid.columns.byName('a17').width = 50;
  // view.grid.columns.byName('a17').cellType = 'html';
  // view.grid.columns.byName('C').width = 50;

  let css = {
    style:{
    textAlign:'center',
    fontWeight:'bold',
    marginTop:'5px'
    }
  };

  let css2 = {
    style:{
    textAlign:'center',
    marginTop:'5px',
    color: "red"
    }
  };

  let css3 = {
    style:{
    textAlign:'center',
    marginTop:'5px',
    color: "green"
    }
  };

  let aaRgx = new RegExp("^[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}$");
  let adRgx = new RegExp("^d[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}$");
  let amRgx = new RegExp("^[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}\\(.+\\)$");
  let admRgx = new RegExp("^d[G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|]{1}\\(.+\\)$$");
  let aNameRgx = new RegExp("^a\\d.*$");
  
  view.grid.onCellPrepare(function (gc) {
    if (gc.isTableCell && gc.gridColumn.name === 'N'){
      gc.style.backColor = 0xff1f77b4;
    }
    if (gc.isTableCell && gc.gridColumn.name === 'C'){
      gc.style.backColor = 0xffff7f00;
    }
    if(gc.isTableCell && !aaRgx.test(gc.cell.value) && !amRgx.test(gc.cell.value) && aNameRgx.test(gc.gridColumn.name)){
      gc.style.element = ui.divText(gc.cell.value, css2);
    }
    if(gc.isTableCell && aaRgx.test(gc.cell.value)){
      console.log(gc)
      gc.style.element = ui.divText(gc.cell.value, css);
    }
    if(gc.isTableCell && amRgx.test(gc.cell.value)){
      console.log(gc)
      gc.style.element = ui.divText(gc.cell.value, css);
    }
    if(gc.isTableCell && adRgx.test(gc.cell.value)){
      console.log(gc)
      gc.style.element = ui.divText(gc.cell.value, css3);
    }
    if(gc.isTableCell && admRgx.test(gc.cell.value)){
      console.log(gc)
      gc.style.element = ui.divText(gc.cell.value, css3);
    }

  });

  pi.close();
}

//name: Peptides
//tags: app
export function Peptides() {

  let appDescription = ui.info(
    [
      ui.span(['For more deatails see LINK ']),
      ui.divText('\n To start the application :', {style: {'font-weight': 'bolder'}}),
      ui.divText('Select the corresponding .csv table with peptide sequences')
    ], 'Transform peptide sequence data to research insights'
  );
  let annotationViewerDiv = ui.div();

  let windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showHelp = false;
  windows.showProperties = false;

  let localTables = grok.shell.tables;
  let namesOfLocalTables = localTables.map((df) => df.name);

  // let chosenFile = ui.choiceInput('File', '', namesOfLocalTables, () => {
  //   main(chosenFile.stringValue, localTables);
  // });

  let chosenFile = ui.button('Open demo', ()=> main());

  let mainDiv = ui.div();
  grok.shell.newView('Peptides', [
    appDescription,
    ui.h2('Choose .csv file'),
    ui.div([
      ui.block25([
        ui.inputs([
          chosenFile
        ])
      ]),
      ui.block75([annotationViewerDiv])
    ]),
    mainDiv
  ]);

  //let table = DG.DataFrame.create();
  //table.name = 'Peptides';
  //let view = grok.shell.addTableView(table);

}