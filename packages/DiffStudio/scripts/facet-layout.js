//name: FacetPlot
//language: javascript
//input: dataframe t {caption: Table}
//input: int colsCount
//input: int max = 4
//input: bool showSplitters = false
//input: bool addIncomplete = true
//input: bool saveGrid = true

const v = grok.shell.addTableView(t);

const cols = t.columns;
const colNames = cols.names();

const node = v.dockManager.dock(DG.Viewer.lineChart(t, {
  xColumnName: colNames[0],
  yColumnNames: colNames.slice(1, colsCount),
  multiAxis: true,
}), 'up', null, 'MultiAxis', 0.7);

const singleLineCharts = [];

const rows = [];
let currentRow = [];

let idx = 0;

for (let i = 1; i < colsCount; ++i) {
  const plot = DG.Viewer.lineChart(t, {
    xColumnName: colNames[0],
    yColumnNames: [colNames[i]],
    autoLayout: false,
    showXAxis: true,
    showYAxis: true,
    showXSelector: false,
    showSplitSelector: false,
    lineWidth: 2,
  });

  ++idx;

  currentRow.push(plot.root);

  if (idx === max) {
    idx = 0;
    const split = ui.splitH(currentRow, null, true);
    rows.push(split);
    currentRow = [];
  }

  singleLineCharts.push(plot);
}

if ((idx > 0) && addIncomplete) {
  if (saveGrid) {
    while (idx < max) {
      ++idx;
      currentRow.push(ui.div(''));
    }
  }

  rows.push(ui.splitH(currentRow, null, true));
}

const facet = ui.splitV(rows, null, true);
const div = ui.div(facet);

grok.shell.dockManager.dock(div, 'fill', node, 'Facet');

facet.style.height = '95%';
facet.style.width = '100%';

const vSplitters = facet.querySelectorAll('div.ui-split-v-divider');
vSplitters.forEach((s) => {s.hidden = !showSplitters;});

const hSplitters = facet.querySelectorAll('div.ui-split-h-divider');
hSplitters.forEach((s) => {s.hidden = !showSplitters;});
