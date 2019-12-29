// Note that the "population" data type becomes int

let t = DataFrame.fromObjects([
    { country: 'USA', population: 321 },
    { country: 'Canada', population: 35 },
    { country: 'Mexico', population: 121 }
]);

grok.addTableView(t);
