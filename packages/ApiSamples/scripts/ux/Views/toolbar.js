//tags: View
//Example of ribbon toolbar, custom header and footer

let header = ui.block([
    ui.divH([
      //ui.iconFA('chart-network'),
      ui.image('https://datagrok.ai/img/logo.svg', 20,20),
      ui.divText('Company name', {style:{margin:'0 5px'}})
    ], {style:{fontWeight:'bolder',alignItems:'center'}})
  ], {style:{padding:'10px', background:'var(--grey-1)'}});
  
  
  let footer = ui.block([
    ui.div([
      'Datagrok example view'
    ], {style:{color:'var(--warm-grey-3)'}})
  ], {style:{padding:'4px', borderTop:'1px solid var(--warm-grey-1)'}});
  
  
  let view = DG.View.create();
  view.name = 'DEMO';
  
  let windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showProperties = false;
  windows.showHelp = false;
  
  let ribbonMenu = DG.Menu.create()
  .group('View settings')
  .item('Top panel', ()=> {
    if (grok.shell.topPanel.innerHTML != '')
      grok.shell.topPanel.innerHTML = ''
    else
      grok.shell.topPanel.append(header)
  })
  .item('Bottom panel', ()=> {
    if (grok.shell.bottomPanel.innerHTML != '')
      grok.shell.bottomPanel.innerHTML = ''
    else
      grok.shell.bottomPanel.append(footer)
  })
  .separator()
  .item('Toolbox', () => windows.showToolbox = !windows.showToolbox)
  .item('Context panel', () => windows.showContextPanel = !windows.showContextPanel)
  .item('Context help', () => windows.showHelp = !windows.showHelp)
  .separator()
  ribbonMenu.item('Hide all', ()=>{
    grok.shell.topPanel.innerHTML = '';
    grok.shell.bottomPanel.innerHTML = '';
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  });
    
  view.ribbonMenu = ribbonMenu; 
  
  let demoData = ui.choiceInput('Data', '', ['', 'Demog', 'Geo', 'Wells'], v=>{
    view.root.innerHTML = '';
    if (v == 'Demog')
      renderStats(demog, v)
      if (v == 'Geo')
        renderStats(geo, v)
        if (v == 'Wells')
          renderStats(wells, v)
  })
  
  view.setRibbonPanels([
    [
      demoData.root
    ],
    [
      ui.tooltip.bind(ui.iconFA('window-maximize', ()=>{
        if (demoData.value == 'Demog')
          grok.shell.newView(demoData.value,[DG.Viewer.grid(demog)])
        if (demoData.value == 'Geo')
          grok.shell.newView(demoData.value,[DG.Viewer.grid(geo)])
        if (demoData.value == 'Wells')
          grok.shell.newView(demoData.value,[DG.Viewer.grid(wells)])
      }), 'Open in new View'),
      ui.tooltip.bind(ui.iconFA('table', ()=>{
        if (demoData.value == 'Demog')
          grok.shell.addTableView(demog)
        if (demoData.value == 'Geo')
          grok.shell.addTableView(geo)
        if (demoData.value == 'Wells')
          grok.shell.addTableView(wells)
      }), 'Open in Table view'),
      ui.tooltip.bind(ui.iconFA('undo', ()=>{
        view.root.innerHTML = '';
        demoData.value = '';
      }), 'Clear view')
    ],
    [
      ui.div([])
    ]
  ])
  
  let demog = grok.data.demo.demog();
  let geo = grok.data.demo.geo();
  let wells = grok.data.demo.wells();
   
  
  function renderStats(data, name){
    view.append(ui.h1(name+' collumns'));
    let stats = ui.divH([],'grok-gallery-grid');
    for (let col of data.columns.toList()) {
      stats.append(ui.div([
        ui.h3(col.name),
        ui.tableFromMap({
          name: col.name,
          totalCount: col.stats.totalCount,
          'missing values': col.stats.missingValueCount,
          valueCount: col.stats.valueCount,
          min: col.stats.min,
          max: col.stats.max,
          avg: col.stats.avg,
          stdev: col.stats.stdev,
          variance: col.stats.variance,
          skew: col.stats.skew,
          kurt: col.stats.kurt,
          med: col.stats.med,
          q1: col.stats.q1,
          q2: col.stats.q2,
          q3: col.stats.q3
        })
      ], 'd4-item-card'))
    }
    view.append(stats);
  }
  
  grok.shell.addView(view);