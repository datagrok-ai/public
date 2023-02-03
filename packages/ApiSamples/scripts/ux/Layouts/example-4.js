//Example of horizontal and verticals splitter layout

let view = grok.shell.newView('Layout 4');
let table = grok.data.testData('demog', 10000);
view.box = true;

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel
windows.showHelp = false; //*hide context help

let tabs = ui.tabControl();

view.append(ui.splitV([
  ui.splitH([
    ui.box(ui.panel([
      ui.inputs([
        ui.columnInput('Column',table),
        ui.stringInput('Process name','Example run 3'),
        ui.choiceInput('Type', '100 ms',['100 ms','200 ms','300 ms']),
        ui.choiceInput('Steps', '10',['10','15','20','30']),
        ui.choiceInput('Control ID', 'AD3343',['AD3343','AD3344','AD3345','AD3346']),
        ui.boolInput('Prediction', false),
        ui.buttonsInput([
          ui.bigButton('Calculate'),
          ui.button('Clear')
        ])
      ]),
    ]),{style:{maxWidth:'400px'}}),
    ui.tabControl({
      'INPUT DATA' : DG.Viewer.fromType('Grid', table).root,
      'SUMMARY' : DG.Viewer.fromType('Statistics', table).root,
    }).root
  ]),
  ui.splitH([
    ui.box(ui.div(DG.Viewer.fromType('Histogram', table, {
      "valueColumnName": "AGE",
      "splitColumnName": "RACE",
      "filteringEnabled": false,
      "showBinSelector": false,
      "showXAxis": true,
      "allowColumnSelection": false,
      "marginLeft": 5,
      "marginRight": 5,
      "xAxisHeight": 15,
      "binWidthRatio": 0.9
	}).root), {style:{maxWidth:'400px'}}),
    DG.Viewer.fromType('Box plot', table).root,
    ui.box(ui.div(DG.Viewer.fromType('Correlation plot', table).root), {style:{maxWidth:'300px'}})
  ], {style:{borderTop:'1px solid var(--grey-2)', marginTop: '5px', paddingTop: '5px'}})
]));