// Different ways to filter, select, and highlight rows

let t = grok.data.demo.demog();
t.rows.filter((row) => row.idx % 3 === 0);
t.rows.select((row) => row.idx % 5 === 0);

let host = ui.divH([
  DG.Viewer.scatterPlot(t, { x: 'height', y: 'weight'}),
  DG.Viewer.histogram(t, { value: 'height' }),
  DG.Viewer.barChart(t),
])

grok.shell.newView('foo').append(host);

let highlightHost = ui.divText('Mouse over to highlight');
highlightHost.onmouseenter = () => t.rows.highlight((i) => i % 2 === 0);
highlightHost.onmouseleave = () => t.rows.highlight(null);

host.appendChild(highlightHost);