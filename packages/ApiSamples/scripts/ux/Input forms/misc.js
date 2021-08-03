//Example app layout with input-form

let view = grok.shell.newView('Input form');
view.box = true;
let windows = grok.shell.windows;
windows.showToolbox = false;
windows.showHelp = false;
windows.showProperties = false;

let query = ui.choiceInput('Query', 'Log Sessions', ['Log Sessions', 'Log Actions', 'Errors Summary On Date', 'Events Summary On Date', 'Unique Users By Date']);
let viewType = ui.choiceInput('Render', 'Scatter plot', ['Scatter plot', 'Bar chart']);
let date = ui.stringInput('Date', 'today');
date.addPatternMenu('datetime');

let filterCol = ui.inputs([
  date,
  query, 
  viewType, 
  ui.buttonsInput([
    ui.bigButton('Run query', ()=>{
      viewerView.innerHTML = '';
      gridView.innerHTML = '';
      grok.shell.info(date.value);
      grok.shell.info('UsageAnalysis:'+query.stringValue.replace(/\s/g, ''));
      grok.data.query('UsageAnalysis:'+query.stringValue.replace(/\s/g, ''), {'date': date.value}).then(t => {
        viewerView.append(DG.Viewer.fromType(viewType.stringValue,t,{'color': 'login'}).root);
        gridView.append(DG.Viewer.grid(t).root);
      });	
    }),
    ui.button('Clear', ()=>{
      viewerView.innerHTML = '';
      gridView.innerHTML = '';
    }),
  ])
]);

let filtersPane = ui.panel([
  ui.h2('Filters'),
  filterCol
]);

let viewerView = ui.box();
let gridView = ui.box();
$(filtersPane).find('input, select').css('min-width','150px');

view.append(ui.splitH([
  ui.box(filtersPane, {style:{maxWidth:'350px'}}),
  ui.splitV([
    viewerView,
    gridView
  ])
]));