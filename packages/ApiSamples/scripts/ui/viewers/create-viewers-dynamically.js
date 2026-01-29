// creating viewers dynamically

let t = grok.data.demo.demog();

// dynamically
let scatterPlot = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t);

// fluent API
let barChart = t.plot.bar();

ui.divV([scatterPlot.root, barChart.root]);