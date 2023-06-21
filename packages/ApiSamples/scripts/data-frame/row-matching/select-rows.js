// Selecting or filtering rows by predicates

let demog = grok.data.demo.demog();
demog.rows.select((row) => row.sex === 'M');
demog.rows.filter((row) => row.age > 42);
grok.shell.addTableView(demog);