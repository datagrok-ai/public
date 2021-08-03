//Example of filtering dataFrame using the Filter Viewer and custom input fields

let view = grok.shell.newView('Filtering');
let table = grok.data.testData('demog');
view.box=true;
let windows = grok.shell.windows;
windows.showToolbox = false;
windows.showHelp = false;
windows.showProperties = false;

let filters = DG.Viewer.fromType('Filters', table, {
  'columnNames':['SEX','AGE','RACE'],
  'showContextMenu':false,
});

let lineChart = DG.Viewer.fromType('Line chart', table, {
  "xColumnName": "AGE",
    "yColumnNames": [
      "WEIGHT"
    ],
    "yAggrTypes": [
      "avg"
    ],
  "splitColumnName": "",
  "segmentColumnName": "SEX",
  "interpolation": "Spline",
  "chartTypes": [
    "Line Chart",
    "Line Chart",
    "Line Chart"
  ],
  "innerChartMarginTop": 0,
  "outerChartMarginRight": 0,
  "outerChartMarginBottom": 0,
  "showXSelector": false,
  "showYSelectors": false,
  "showAggrSelectors": false,
  "showCloseLink": false,
  "showSplitSelector": false,
  "showYAxis": false,
  "legendVisibility": "Never",
  "showMarkers": "Never",
});

let barChart = DG.Viewer.fromType('Bar chart', table, {
  "splitColumnName": "Site",
  "stackColumnName": "Sex",
  "valueColumnName": "AGE",
  "valueAggrType": "count",
  "legendVisibility": "Never",
  "barSortType": "by value",
  "barSortOrder": "desc",
  "showValueAxis": false,
  "showValueSelector": false,
  "showCategorySelector": false,
  "showStackSelector": false,
  "showSelectedRows": false,
  'Title':'Bar chart'
});

barChart.root.style.height = '150px';
lineChart.root.style.height = '150px';
filters.root.style.overflow = 'scroll';

view.append(ui.splitH([
ui.splitV([
  // Add custom filter  
  ui.panel([
    ui.h1('Custom filters'),
    ui.inputs([
      ui.choiceInput('Sex', '',['M', 'F'], (v)=>{
        table.rows.match('Sex ='+v).filter();
      }),
      ui.choiceInput('Race', '',['Caucasian', 'Asian', 'Black', 'Other'], (v)=>{
        table.rows.match('Race = '+v).filter();
      }),
      ui.intInput('Age', '', (v)=>{
        table.rows.match('Age ='+v).filter();
      }),
      ui.stringInput('Weight', '', (v)=>{
        table.rows.match('Weight >'+v).filter();
      }),
    ])
  ]),
  // Add filter viewer
  ui.box(ui.panel([ui.h1('Default filters')]), {style:{maxHeight:'55px'}}),
  filters.root
], {style:{maxWidth:'270px'}}),
ui.splitV([
  ui.box(
  ui.panel([
    ui.block(ui.h1('Demography insights')),
    ui.block50(barChart.root),
    ui.block50(lineChart.root),
  ]), {style:{maxHeight:'220px'}}),
  ui.box(ui.panel(ui.h1('Demography data')), {style:{maxHeight:'55px'}}),
  DG.Viewer.fromType('Grid',table).root
])
]));