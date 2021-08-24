//tags: View
//Toolbox

let view = grok.shell.newView('Toolbox simple view');
view.box = true;

let rows = ui.intInput('Rows', 100);
let choiceData = ui.choiceInput('Data', '',['','demog','wells','geo']);
let data = '';
let chartView = ui.splitV([]); 
let gridView = ui.box();
let stats = ui.div();
let statsView = ui.box();
statsView.style.maxWidth='180px';
statsView.style.maxHeight = '280px';
let colInput;

view.append(ui.splitH([gridView, chartView]))
let scatter = ui.iconFA('chart-scatter',()=>{
  if (data != '')
  	chartView.append(DG.Viewer.fromType('Scatter plot', data).root)
});
scatter.style.fontSize = '20px';
scatter.style.fontWeight = '500';
scatter.style.color = 'var(--warm-grey-3)';

let lineChart = ui.iconFA('chart-line',()=>{
  if (data != '')
  	chartView.append(DG.Viewer.fromType('Line chart', data).root)
});
lineChart.style.fontSize = '20px';
lineChart.style.fontWeight = '500';
lineChart.style.color = 'var(--warm-grey-3)';

let barChart = ui.iconFA('chart-bar',()=>{
  if (data != '')
  	chartView.append(DG.Viewer.fromType('Bar chart', data).root)
});
barChart.style.fontSize = '20px';
barChart.style.fontWeight = '500';
barChart.style.color = 'var(--warm-grey-3)';

let pieChart = ui.iconFA('chart-pie',()=>{
  if (data != '')
  	chartView.append(DG.Viewer.fromType('Pie chart', data).root)
});
pieChart.style.fontSize = '20px';
pieChart.style.fontWeight = '500';
pieChart.style.color = 'var(--warm-grey-3)';

let acc = ui.accordion();
acc.addPane('Demo data',()=>ui.buttonsInput([
    choiceData,
    rows,
      ui.bigButton('Apply',()=>{
		if (choiceData.value == 'demog')
		  data = grok.data.demo.demog(rows.value)
		if (choiceData.value == 'wells')
		  data = grok.data.demo.wells(rows.value)
		if (choiceData.value == 'geo')
		  data= grok.data.demo.geo(rows.value)
		if (choiceData.value == ''){
		  data = '';
          gridView.innerHTML = '';
          chartView.innerHTML = '';
          stats.innerHTML = '';
		  grok.shell.error('Select data');
        }  
        gridView.innerHTML = '';
        chartView.innerHTML = '';
        if (data != ''){
        colInput = ui.choiceInput('Column', '', data.columns.toList(), col=>{
          statsView.innerHTML = '';
          statsView.append(ui.tableFromMap({
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
          }))
        });
        gridView.append(DG.Viewer.grid(data).root);
        stats.innerHTML = '';
        statsView.innerHTML = '';
        stats.append(ui.divV([colInput,statsView]));
        }
      })
]), true);

acc.addPane('Viewers',()=>ui.divV([
  ui.divH([
    ui.tooltip.bind(scatter, 'Scatter plot'),
    ui.tooltip.bind(lineChart, 'Line chart'),
    ui.tooltip.bind(barChart, 'Bar chart'),
    ui.tooltip.bind(pieChart, 'Pie chart'),
  ], {style:{marginRight:'20px',justifyContent:'space-between'}})
]));
acc.addPane('Statistic',()=> stats);
acc.addPane('Actions',()=>ui.divV([
  ui.link('Add View', ()=>{grok.shell.addTableView(data)}, '',''),
  ui.link('Add Filters',()=>{grok.shell.dockManager.dock(DG.Viewer.fromType('Filters', data).root, 'left',null,'Filters',0.2);},'',''),
  ui.link('General Info', ()=>{
    grok.shell.info(data.name);
    grok.shell.info(data.rowCount);
  }, '',''),
  ui.link('Clear View', ()=>{
    data = '';
    gridView.innerHTML = '';
    chartView.innerHTML = '';
    stats.innerHTML = '';
  }, '','')
]));

view.toolbox = acc.root