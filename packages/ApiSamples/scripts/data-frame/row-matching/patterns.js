// Matching rows by patterns.
//

let d = grok.data.demo.demog(20000);
grok.shell.addTableView(d).scatterPlot();

d.rows.match({ age: '>30', sex: 'M'}).filter();
d.rows.match('age > 50').select();
d.rows.match('race = Asian').highlight();