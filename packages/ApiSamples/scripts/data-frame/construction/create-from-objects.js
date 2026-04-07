// Note that the "population" data type becomes int

let t = DG.DataFrame.fromObjects([
  {country: 'USA', population: 321},
  {country: 'Canada', population: 35},
  {country: 'Mexico', population: 121}
]);

grok.shell.addTableView(t);
