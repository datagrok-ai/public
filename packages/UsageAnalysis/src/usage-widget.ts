/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";

export class UsageWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    super(ui.box());
   
    let Summary: any = {
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

    grok.data
      .query('UsageAnalysis:NewUsersEventsErrors')
      .then((t) => {
        for (let i = 0; i < 3; i++) {
          if (t.get('counter_type', i) == 'users_count') {
            Summary.users.today = t.get('day', i);
            Summary.users.week = t.get('week', i);
            Summary.users.month = t.get('month', i)
          }
          if (t.get('counter_type', i) == 'events_count') {
            Summary.events.today = t.get('day', i);
            Summary.events.week = t.get('week', i);
            Summary.events.month = t.get('month', i)
          }
          if (t.get('counter_type', i) == 'errors_count') {
            Summary.errors.today = t.get('day', i);
            Summary.errors.week = t.get('week', i);
            Summary.errors.month = t.get('month', i)
          }
        }
        userStats.appendChild(Summary.users.month);
        eventsStats.appendChild(Summary.events.month);
        errorsStats.appendChild(Summary.errors.month);
      });

    let period = ui.choiceInput('', 'Month', ['Today','Week','Month'], (v: string) => {
      const charts = [usersChart, eventsChart, errorsChart, userStats, eventsStats, errorsStats];
      let days: Record<string, string> = {
        'Today': '2',
        'Week': '7',
        'Month': '31'
      }

      charts.forEach((c) => c.innerHTML = '');
      usersChart.appendChild(getNewUsers(days[v]));
      eventsChart.appendChild(getNewEvents(days[v]));
      errorsChart.appendChild(getNewErrors(days[v]));
      userStats.appendChild(Summary.users[v.toLowerCase()])
      eventsStats.appendChild(Summary.events[v.toLowerCase()])
      errorsStats.appendChild(Summary.errors[v.toLowerCase()])
    });

    let usersChart = ui.box();
    usersChart.style.margin = '0px 12px';
    let userStats = ui.div([], {style: {color: 'var(--grey-4)'}});
    let users = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Users', {style: {marginBottom: '8px', fontWeight: 'bold'}}),
        userStats
      ]), {style: {maxWidth: '80px'}}),
      usersChart
    ]);

    let eventsChart = ui.box();
    eventsChart.style.margin = '0px 12px';
    let eventsStats = ui.div([], {style: {color: 'var(--grey-4)'}});
    let events = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Events', {style: {marginBottom: '8px', fontWeight: 'bold'}}),
        eventsStats
      ]), {style: {maxWidth: '80px'}}),
      eventsChart
    ]);

    let errorsChart = ui.box();
    errorsChart.style.margin = '0px 12px';
    let errorsStats = ui.div([], {style: {color: 'var(--grey-4)'}});
    let errors = ui.splitH([
      ui.box(ui.panel([
        ui.divText('Errors', {style: {marginBottom: '8px', fontWeight: 'bold'}}),
        errorsStats
      ]), {style: {maxWidth: '80px'}}),
      errorsChart
    ]);

  
    usersChart.appendChild(getNewUsers('30'));
    eventsChart.appendChild(getNewEvents('30'));
    errorsChart.appendChild(getNewErrors('30'));

    this.root.appendChild(ui.splitV([
      ui.box(
        ui.divH([
          period.root,
          ui.link('Open Usage Analysis',()=>{grok.functions.call("UsageAnalysis:usageAnalysisApp")},'Open Usage Analysis application')        
        ], {style:{justifyContent:'space-between', alignItems:'center', margin:'0 12px'}})
      , {style:{maxHeight:'40px', marginBottom:'12px'}}),
      users,
      events,
      errors,
    ]));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Usage');
    this.order = super.addProperty('order', DG.TYPE.STRING, '2');
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
          root.appendChild(chart1.root)
        });
  return root;      
}

function getNewEvents(days:string){
  let root = ui.box();
  grok.data.query('UsageAnalysis:NewEventsSummaryByDate', {'days':days})
    .then((t) => {
        let chart2 = DG.Viewer.fromType('Line chart', t, eventsChartStyle);
        chart2.root.style.maxHeight = '85px';
        root.appendChild(chart2.root)
      }); 
  return root;      
}

function getNewErrors(days:string){
  let root = ui.box();
  grok.data.query('UsageAnalysis:NewErrorsSummaryByDate', {'days':days})
    .then((t) => {
        let chart3 = DG.Viewer.fromType('Line chart', t, errorsChartStyle);
        chart3.root.style.maxHeight = '85px';
        root.appendChild(chart3.root)
      });
  return root;      
}

