//tags: View
//Context panel

let testData = {
    demog: grok.data.testData('demog'),
    wells: grok.data.testData('wells'),
    biosensor: grok.data.testData('biosensor'),
    walk: grok.data.testData('random walk'),
    geo: grok.data.testData('geo'),
  }
  
  let viewerOptions = {
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
    };
  
  class Demo {
    constructor(name, description, dataFrame) {
      this.name = name;
      this.description = description;
      this.dataFrame = dataFrame;
    }
  }
  
  class DemoHandler extends DG.ObjectHandler {
    get type() { return 'demo' }
    isApplicable(x) { return x instanceof Demo; }
    
    renderProperties(x) { 
      let chart = DG.Viewer.fromType('Histogram', x.dataFrame, viewerOptions);
      let grid = DG.Viewer.fromType('Grid', x.dataFrame);
      chart.root.style.maxWidth = '300px';
      grid.root.style.maxWidth = '300px';
      
      let acc = ui.accordion();
      acc.addPane('Details', ()=>
        ui.tableFromMap({
          Table: x.dataFrame.name,
          Rows: x.dataFrame.rowCount,
          Columns: x.dataFrame.columns.length
        })
      );
      acc.addPane('Preview', ()=> ui.div(grid.root));
      acc.addPane('Code', ()=> ui.divText(`let data = grok.data.testData('${x.name}', 1000); \n grok.shell.addTableView(data);`));
      acc.addPane('Actions', ()=>ui.divV([
        ui.link('Add to view',()=>{
          grok.shell.dockManager.dock(DG.Viewer.fromType('Grid', x.dataFrame).root, 'down', grok.shell.v.dockNode ,x.name, 0.5)
        },'Dock grid to current view', ''),
        ui.link('Open table', ()=>{
          grok.shell.addTableView(x.dataFrame);
        }, 'Open table and table view'),
        ui.link('Delete', ()=>{
         $('body').find('.item-'+`${x.name}`).parent().detach();
        }, 'Delete from current view', '')
      ]));
      
      return ui.panel([
        ui.h1(`${x.name}`),
        ui.divText(`${x.description}`),
        chart.root,
        acc.root
      ]); 
    }
    renderTooltip(x) { return ui.divText(`${x.description}, \n Click to see object details`); }
    
    renderCard(x) {
      return ui.bind(x, ui.divV([
        ui.h2(`${x.name}`),
        ui.divText(`${x.description}`),
          ui.tableFromMap({
          'Rows':x.dataFrame.rowCount,
        })
      ], `item-${x.name}`));
    }
    
  }
  
  DG.ObjectHandler.register(new DemoHandler());
  let demog = new Demo('Demography', 'Clinical study demographics data', testData.demog);
  let wells = new Demo('Wells', 'Experimental plate wells: barcode, row, col, pos, etc', testData.wells);
  let biosensor = new Demo('Biosensor', 'Wearable sensor data: time, x, y, z, temp, eda', testData.biosensor);
  let walk = new Demo('Random-walk', 'Random walk data for the specified number of dimensions ', testData.walk);
  let geo = new Demo('Geo', 'Geographic coordinates given as latitude/longitude pairs', testData.geo);
  
  
  let windows = grok.shell.windows;
  windows.showHelp = false;
  
  let view = grok.shell.newView('Context panel');
  view.append(ui.divV([
    ui.panel([
      ui.h1('Datagrok test data'),
      ui.divText('Generic dataset with the defined number of rows and columns'),
    ]),
    ui.divH([
      ui.renderCard(demog),
      ui.renderCard(wells),
      ui.renderCard(biosensor),
      ui.renderCard(walk),
      ui.renderCard(geo)
    ], 'grok-gallery-grid')
  ]));
  