// Matching rows by patterns.

let d = grok.data.demo.demog(1000);
grok.shell.addTableView(d).scatterPlot();

d.rows.match('age > 50').select();
d.rows.match('race = Asian').highlight();