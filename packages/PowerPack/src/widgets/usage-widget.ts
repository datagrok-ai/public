/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";
import {divText} from "datagrok-api/ui";

export class UsageWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.box());
   
    let Summary = {
      users: {
        today:'',
        week:'',
        month:''
      },
      events: {
        today:'',
        week:'',
        month:''
      },
      errors: {
        today:'',
        week:'',
        month:''
      }
    };

    grok.data.query('UsageAnalysis:NewUsersEventsErrors')
    .then((t) => {
      for (let i = 0; i<3; i++){
        if (t.get('counter_type', i) == 'users_count'){
          Summary.users.today = t.get('day', i);
          Summary.users.week = t.get('week', i);
          Summary.users.month = t.get('month', i)
        }
        if (t.get('counter_type', i) == 'events_count'){
          Summary.events.today = t.get('day', i);
          Summary.events.week = t.get('week', i);
          Summary.events.month = t.get('month', i)
        }
        if (t.get('counter_type', i) == 'errors_count'){
          Summary.errors.today = t.get('day', i);
          Summary.errors.week = t.get('week', i);
          Summary.errors.month = t.get('month', i)
        }
      }
      userStats.append(Summary.users.month);
      eventsStats.append(Summary.events.month);
      errorsStats.append(Summary.errors.month);
    })

    let period = ui.choiceInput('', 'Month', ['Today','Week','Month'], (v:string) => {
      switch (v){
        case 'Today':
          console.log('Today');
          usersChart.innerHTML = '';
          usersChart.append(getNewUsers('2'));
          eventsChart.innerHTML = '';
          eventsChart.append(getNewEvents('2'));
          errorsChart.innerHTML = '';
          errorsChart.append(getNewErrors('2'));
          userStats.innerHTML = '';
          userStats.append(Summary.users.today)
          eventsStats.innerHTML = '';
          eventsStats.append(Summary.events.today)
          errorsStats.innerHTML = '';
          errorsStats.append(Summary.errors.today)
          break;
        case 'Week':
          console.log('Week');
          usersChart.innerHTML = '';
          usersChart.append(getNewUsers('7'));
          eventsChart.innerHTML = '';
          eventsChart.append(getNewEvents('7'));
          errorsChart.innerHTML = '';
          errorsChart.append(getNewErrors('7'));
          userStats.innerHTML = '';
          userStats.append(Summary.users.week)
          eventsStats.innerHTML = '';
          eventsStats.append(Summary.events.week)
          errorsStats.innerHTML = '';
          errorsStats.append(Summary.errors.week)
          break;
        case 'Month':
          console.log('Month');
          usersChart.innerHTML = '';
          usersChart.append(getNewUsers('30'));
          eventsChart.innerHTML = '';
          eventsChart.append(getNewEvents('30'));
          errorsChart.innerHTML = '';
          errorsChart.append(getNewErrors('30'));
          userStats.innerHTML = '';
          userStats.append(Summary.users.month)
          eventsStats.innerHTML = '';
          eventsStats.append(Summary.events.month)
          errorsStats.innerHTML = '';
          errorsStats.append(Summary.errors.month)
          break;
      }
        
    })

    let usersChart = ui.box();
    usersChart.style.margin = '0px 12px';
    let userStats = ui.div([],{style:{color:'var(--grey-4)'}});
    let users = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Users', {style:{marginBottom:'8px', fontWeight:'bold'}}),
        userStats
      ]), {style:{maxWidth:'80px'}}),
        usersChart
    ]);

    let eventsChart = ui.box();
    eventsChart.style.margin = '0px 12px';
    let eventsStats = ui.div([],{style:{color:'var(--grey-4)'}});
    let events = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Events', {style:{marginBottom:'8px', fontWeight:'bold'}}),
        eventsStats
      ]), {style:{maxWidth:'80px'}}),
        eventsChart
    ]);

    let errorsChart = ui.box();
    errorsChart.style.margin = '0px 12px';
    let errorsStats = ui.div([],{style:{color:'var(--grey-4)'}});
    let errors = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Errors', {style:{marginBottom:'8px', fontWeight:'bold'}}),
        errorsStats
      ]), {style:{maxWidth:'80px'}}),
        errorsChart
    ]);

  
    usersChart.append(getNewUsers('30'));
    eventsChart.append(getNewEvents('30'));
    errorsChart.append(getNewErrors('30'));

    this.root.append(ui.splitV([
      ui.box(
        ui.divH([
          period.root,
          ui.link('Open Usage Analysis',()=>{grok.functions.call("UsageAnalysis:startApp")},'Open Usage Analysis application')        
        ], {style:{justifyContent:'space-between', alignItems:'center', margin:'0 12px'}})
      , {style:{maxHeight:'40px', marginBottom:'12px'}}),
      users,
      events,
      errors,
    ]));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Usage');
  }
}

let usersChartStyle = {
  "aggrType": "count",
  "innerChartMarginTop": 0,
  "innerChartMarginBottom": 0,
  "outerChartMarginTop": 5,
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
  "Title":"Users"
};

let eventsChartStyle = {
  "defaultLineColor": 4281449313,
  "innerChartMarginTop": 0,
  "innerChartMarginBottom": 0,
  "outerChartMarginTop": 5,
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
  "innerChartMarginBottom": 0,
  "outerChartMarginTop": 5,
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

function getNewUsers(days:string){
  let root = ui.box();
  grok.data.query('UsageAnalysis:NewUsersSummaryByDate', {'days':days})
        .then((t) => {
          let chart1 = DG.Viewer.fromType('Line chart', t, usersChartStyle);
          chart1.root.style.maxHeight = '85px';
          root.append(chart1.root)
        });
  return root;      
}

function getNewEvents(days:string){
  let root = ui.box();
  grok.data.query('UsageAnalysis:NewEventsSummaryByDate', {'days':days})
    .then((t) => {
        let chart2 = DG.Viewer.fromType('Line chart', t, eventsChartStyle);
        chart2.root.style.maxHeight = '85px';
        root.append(chart2.root)
      }); 
  return root;      
}

function getNewErrors(days:string){
  let root = ui.box();
  grok.data.query('UsageAnalysis:NewErrorsSummaryByDate', {'days':days})
    .then((t) => {
        let chart3 = DG.Viewer.fromType('Line chart', t, errorsChartStyle);
        chart3.root.style.maxHeight = '85px';
        root.append(chart3.root)
      });
  return root;      
}

