/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";
import {divText} from "datagrok-api/ui";

export class UsageWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.box());
    
    let eventsChartStyle = {
      "defaultLineColor": 4281449313,
      "innerChartMarginTop": 0,
      "outerChartMarginRight": 0,
      "outerChartMarginBottom": 0,
      "yGlobalScale": false,
      "showTopPanel": false,
      "showMouseOverRowLine": false,
      "showXSelector": false,
      "showYSelectors": false,
      "showAggrSelectors": false,
      "showSplitSelector": false,
      "showYAxis": false,
      "showMarkers": "Never",
    };

    

    let usersChartStyle = {
      "aggrType": "count",
      "innerChartMarginTop": 0,
      "outerChartMarginRight": 0,
      "outerChartMarginBottom": 0,
      "yGlobalScale": false,
      "showTopPanel": false,
      "showMouseOverRowLine": false,
      "showXSelector": false,
      "showYSelectors": false,
      "showAggrSelectors": false,
      "showSplitSelector": false,
      "showYAxis": false,
      "showMarkers": "Never",
    };

    let errorsChartStyle = {
      "defaultLineColor": 4294262794,
      "innerChartMarginTop": 0,
      "outerChartMarginRight": 0,
      "outerChartMarginBottom": 0,
      "yGlobalScale": false,
      "showTopPanel": false,
      "showMouseOverRowLine": false,
      "showXSelector": false,
      "showYSelectors": false,
      "showAggrSelectors": false,
      "showSplitSelector": false,
      "showYAxis": false,
      "showMarkers": "Never",
    };

    grok.data.query('UsageAnalysis:NewUsersEventsErrors')
    .then((t) => {
      for (let i = 0; i<3; i++){
        if (t.get('counter_type', i) == 'users_count')
          userStats.append(ui.tableFromMap({
            Today: t.get('day', i),
            Month: t.get('month', i),
          }))
        if (t.get('counter_type', i) == 'events_count')
          eventsStats.append(ui.tableFromMap({
            Today: t.get('day', i),
            Month: t.get('month', i),
          })) 
        if (t.get('counter_type', i) == 'errors_count')
          errorsStats.append(ui.tableFromMap({
            Today: t.get('day', i),
            Month: t.get('month', i),
          })) 
      }
    })

    let usersChart = ui.box();
    usersChart.style.marginLeft = '25px';
    let userStats = ui.div();
    let users = ui.splitH([
      ui.box(ui.div([
        ui.divText('Users', {style:{fontWeight:'bold', marginLeft:'12px',marginBottom:'-4px'}}),
        userStats
      ]), {style:{maxWidth:'120px'}}),
        usersChart
    ]);

    let eventsChart = ui.box();
    eventsChart.style.marginLeft = '25px';
    let eventsStats = ui.div();
    let events = ui.splitH([
      ui.box(ui.div([
        ui.divText('Events', {style:{fontWeight:'bold',marginLeft:'12px',marginBottom:'-4px'}}),
        eventsStats
      ]), {style:{maxWidth:'120px'}}),
        eventsChart
    ]);

    let errorsChart = ui.box();
    errorsChart.style.marginLeft = '25px';
    let errorsStats = ui.div();
    let errors = ui.splitH([
      ui.box(ui.div([
        ui.divText('Errors', {style:{fontWeight:'bold', marginLeft:'12px',marginBottom:'-4px'}}),
        errorsStats
      ]), {style:{maxWidth:'120px'}}),
        errorsChart
    ]);

    grok.data.query('UsageAnalysis:NewUsersSummaryByDate', {'days':'30'})
        .then((t) => {
          let chart1 = DG.Viewer.fromType('Line chart', t, usersChartStyle);
          chart1.root.style.maxHeight = '80px';
          usersChart.append(chart1.root)
        });

    grok.data.query('UsageAnalysis:NewEventsSummaryByDate', {'days':'30'})
    .then((t) => {
        let chart2 = DG.Viewer.fromType('Line chart', t, eventsChartStyle);
        chart2.root.style.maxHeight = '80px';
        eventsChart.append(chart2.root)
      }); 
      
    grok.data.query('UsageAnalysis:NewErrorsSummaryByDate', {'days':'30'})
    .then((t) => {
        let chart3 = DG.Viewer.fromType('Line chart', t, errorsChartStyle);
        chart3.root.style.maxHeight = '80px';
        errorsChart.append(chart3.root)
      });    
    
    this.root.append(ui.splitV([
      users,
      events,
      errors,
      ui.box(
        ui.div([ui.link('Open Usage Analysis',()=>{grok.functions.call("UsageAnalysis:startApp")},'Open Usage Analysis application')        
        ]), {style:{maxHeight:'20px', marginLeft:'12px'}})
    ], {style:{marginTop:'10px'}} ));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Usage');
  }
}

