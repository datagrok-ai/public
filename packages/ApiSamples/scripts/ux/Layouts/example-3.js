//Splitter with Tabs

let view = grok.shell.newView('Layout 3');
let table = grok.data.testData('demog', 10000);
view.box = true;

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel
windows.showHelp = false; //*hide context help

let sideChart = ui.block();
let sideDetails = ui.block();
let sideAcc = ui.accordion();

let covidTop = ui.panel();
let covidTab = ui.splitV([covidTop]);

let crimeTab = ui.box();

let dataHealth = grok.data.getDemoTable("geo/health-and-wealth.csv").then((t)=>{
  let chart = DG.Viewer.fromType('Histogram', t);
  chart.setOptions({
    "valueColumnName": "PopulationMarkerSize",
    "splitColumnName": "Region",
    "bins": 33,
    "filteringEnabled": false,
    "showCurrentRow": false,
    "showMouseOverRow": false,
    "showFilteredOutRows": false,
    "showBinSelector": false,
    "showColumnSelector": false,
    "showRangeSlider": false,
    "showMouseOverRowGroup": false,
    "showXAxis": true,
    "showContextMenu": false,
    "allowColumnSelection": false,
    "rowIndicatorSize": 32,
    "marginLeft": 5,
    "marginRight": 5,
    "filterMarginBottom": 27,
    "binWidthRatio": 0.9,
    "xAxisHeight": 30,
    "allowDynamicMenus": false,
    "title": "",
    "description": " ",
    "descriptionPosition": "Center",
    "descriptionVisibilityMode": "Always"  
  }); 
  chart.root.style.maxHeight = '150px';
  sideChart.append(ui.divText('Popilation marker size'));
  sideChart.append(chart.root);
  
  let chart2 = DG.Viewer.fromType('Histogram', t);
  chart2.setOptions({
    "valueColumnName": "Life Expectancy At Birth",
  "colorColumnName": "PopulationMarkerSize",
  "filteringEnabled": false,
  "showCurrentRow": false,
  "showMouseOverRow": false,
  "showFilteredOutRows": false,
  "showBinSelector": false,
  "showColumnSelector": false,
  "showRangeSlider": false,
  "showMouseOverRowGroup": false,
  "showXAxis": true,
  "showContextMenu": false,
  "allowColumnSelection": false,
  "marginLeft": 5,
  "marginRight": 5,
  "binWidthRatio": 0.9,
  "allowDynamicMenus": false,
  "description": "Life expectancy at birth",
  "descriptionPosition": "Center",
  "descriptionVisibilityMode": "Always"  
  });
  chart2.root.style.maxHeight = '150px';
  sideChart.append(chart2.root);
  sideAcc.addPane('Details',()=>
    ui.tableFromMap({
      Table: t.name,
      Rows: t.rowCount,
      Columns: t.columns.length
  	})
  )
  sideAcc.addPane('Countries details',()=>
    ui.tableFromMap({
      Countires: t.columns.byName('Country').stats.uniqueCount,
      'Values count': t.columns.byName('Country').stats.totalCount,
      'Type': t.columns.byName('Country').type,
      'First value': t.columns.byName('Country').categories[0],
      'Last value': t.columns.byName('Country').categories[251],
  	})
  );
  sideAcc.addPane('Country list',()=>
    ui.list(t.columns.byName('Country').categories)
  );
  sideDetails.append(ui.h1('Data information'));
  sideDetails.append(sideAcc.root); 
  sideDetails.append(ui.panel([
    ui.bigButton('Export data',()=> grok.shell.info('Success!')),
    ui.button(ui.iconFA('redo'),()=> grok.shell.info('Success!'))
  ])); 
});

let dataCovid = grok.data.getDemoTable("covid-19-cases.csv").then((t)=>{
  let chart = DG.Viewer.fromType('Histogram', t);
  chart.setOptions({
  "valueColumnName": "Date",
  "bins": 40,
  "showCurrentRow": false,
  "showBinSelector": false,
  "showColumnSelector": false,
  "showRangeSlider": false,
  "showXAxis": true,
  "showYAxis": false,
  "showContextMenu": false,
  "allowColumnSelection": false,
  "marginLeft": 5,
  "marginRight": 5,
  "binWidthRatio": 0.9,
  "allowDynamicMenus": false
  });
  chart.root.style.height = '200px';
  covidTop.append(ui.block([
    ui.block25([
      ui.label('Total'),
      ui.divText(t.columns.byName('Value').stats.sum.toString(), {style:{fontSize:'25px', fontWeight:'bold',margin:'5px 0'}})
    ]),
    ui.block25([
      ui.label('Count'),
      ui.divText(t.columns.byName('Value').stats.valueCount.toString(), {style:{fontSize:'25px', fontWeight:'bold',margin:'5px 0'}})
    ]),
    ui.block25([
      ui.label('Maximum'),
      ui.divText(t.columns.byName('Value').max.toString(), {style:{fontSize:'25px', fontWeight:'bold',margin:'5px 0'}})
    ]),
    ui.block25([
      ui.label('Minimun'),
      ui.divText(t.columns.byName('Value').min.toString(), {style:{fontSize:'25px', fontWeight:'bold',margin:'5px 0'}})
    ])
  ], {style:{marginBottom:'20px'}}));
  covidTop.append(ui.block(chart.root));
  covidTab.append(ui.block(DG.Viewer.fromType('Grid', t).root));
});

let accDemog = ui.accordion();
accDemog.addPane('Scatter plot',()=> ui.block([
  ui.h1('Scatter plot'),
  DG.Viewer.fromType('Scatter plot', table, {}).root
]), true);
accDemog.addPane('Bar charts',()=> ui.divH([
  ui.block([
    ui.h1('Bar chart 1'),
    DG.Viewer.fromType('Bar chart', table, {'splitColumnName':'Race'}).root
  ]),
  ui.block([
    ui.h1('Bar chart 2'),
    DG.Viewer.fromType('Bar chart', table, {'splitColumnName':'Sex'}).root
  ]),
  ui.block([
    ui.h1('Bar chart 3'),
    DG.Viewer.fromType('Bar chart', table, {}).root
  ])
]));
accDemog.addPane('Data grid',()=> ui.block([
  ui.h1('Data grid'),
  DG.Viewer.fromType('Grid', table, {}).root
]));
accDemog.addPane('Statistics',()=> ui.block([
  ui.h1('Statistics'),
  DG.Viewer.fromType('Statistics', table, {}).root
]));

let dataCrime = grok.data.getDemoTable("crime.csv").then((t)=>{
  crimeTab.append(DG.Viewer.fromType('Grid', t).root);
});
      
view.append(ui.splitH([
	ui.box(
      ui.panel([
        ui.h1('Summary'),  
        sideChart,
        sideDetails,
    ]),{style:{maxWidth:'25%'}}),
  	ui.tabControl({
      'COVID 19':()=>covidTab,
      'DEMOG':()=> ui.panel(accDemog.root),
      'CRIME':()=>crimeTab,
    }).root
])
);
