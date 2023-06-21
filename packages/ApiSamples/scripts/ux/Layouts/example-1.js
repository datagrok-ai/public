//Splitter layout

let data = grok.data.demo.demog();
let view = grok.shell.newView('Layout 1');
view.box = true;

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel
windows.showHelp = false; //*hide context help

//inputs
let colInput = ui.columnInput('Column', data, data.col('age'), ()=>{
});

let searchInput = ui.searchInput('', '', ()=>{
  if (searchInput.value != '')
    data.rows.match(colInput.stringValue+'='+searchInput.value).filter();
});
let filters = ui.button(ui.iconFA('filter'),()=>{
  grok.shell.dockManager.dock(filterDock.root, 'left',null,'Filters',0.2);
});

//viewers
let barChart = DG.Viewer.barChart(data);
barChart.setOptions({
  split:'race',
});
let boxPlot = DG.Viewer.boxPlot(data);
boxPlot.setOptions({
  split:'weight'
});
let grid = DG.Viewer.grid(data);
let statistics = DG.Viewer.fromType('Statistics', data);
let filterDock = DG.Viewer.fromType('Filters', data);

//layout
view.append(ui.splitV([
  ui.box(ui.panel([
    ui.divH([colInput,searchInput,filters])
  ]), {style:{maxHeight:'80px'}}),
  ui.splitH([
    barChart.root,
    boxPlot.root,
  ]),
  ui.tabControl({
    'DATA': ()=> grid.root,
    'STATISTICS': ()=> statistics.root
  }).root
]));