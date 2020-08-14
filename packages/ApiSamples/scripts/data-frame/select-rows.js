// Selecting or filtering rows by predicates

let demog = grok.data.testData('demog', 5000);
demog.rows.select((row) => row.sex === 'M');
demog.rows.filter((row) => row.age > 42);
grok.shell.add(demog);