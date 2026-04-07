// Creating HTML-driven tables
// For high-performance grid, check out samples under /grid

let users = [
  {name: 'Mike', surname: 'Moore', age: 30, sex: 'Male'},
  {name: 'Kate', surname: 'Jackson', age: 27, sex: 'Female'},
  {name: 'John', surname: 'Stone', age: 32, sex: 'Male'},
];

let table1 = ui.tableFromMap({
  user: grok.shell.user.toMarkup(),
  project: grok.shell.project.toMarkup(),
  time: new Date(),
});

let table2 = DG.HtmlTable.create(users, (item) => [item.name, item.surname, item.age, item.sex]);

let table3 = ui.table(users, (item) => [item.name, item.surname, item.age, item.sex], 
  ['First name', 'Last name', 'Age', 'Sex']);

let view = grok.shell.newView('table with renderer', [
  table1,
  table2,
  table3,
]);

table2.remove(users[2]);