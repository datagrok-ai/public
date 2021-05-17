// Different ways to get a column

let d = grok.data.demo.demog();

let ageColumn = d.getCol('age');
ageColumn = d.columns.byName('age');
ageColumn = d.columns['age'];
ageColumn = d.columns.byIndex(3);
ageColumn = d.columns[3];