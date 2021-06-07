// Creating HTML-driven tables
// For high-performance grid, check out samples under /grid

grok.shell.newView('table from map', [
  ui.tableFromMap({
    user: grok.shell.user.toMarkup(),
    project: grok.shell.project.toMarkup(),
    time: new Date(),
  })
]);

let view = grok.shell.newView('table with renderer');
var myList = [
  {key: 'first object', value: 5},
  {key: 'second object', value: false},
  {key: 'third object', value: false},
];

let table = DG.HtmlTable.create(myList, (item, idx) => [item.key, item.value]);

//and let's remove the third object from the table
table.remove(myList[2]);

view.root.appendChild(table.root);