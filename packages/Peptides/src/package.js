/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();
var alignementModule = null;

function main(tableName, localTables) {
  grok.shell.info(tableName);

  let peptides = localTables.find(({name}) => name === tableName);
  
  let view = grok.shell.addTableView(peptides);

  grok.shell.v = view;

}


//name: Peptides
//tags: app
export function Peptides() {

  alignementModule = await Module();

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
  let chosenFile = ui.choiceInput('File', '', namesOfLocalTables, () => {
    main(chosenFile.stringValue, localTables);
  });

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