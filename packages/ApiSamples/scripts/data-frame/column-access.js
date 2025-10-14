//tags: Column, access
// Different ways to get a column

let d = grok.data.demo.demog();

let ageColumn = d.getCol('age');
ageColumn = d.columns.byName('age');
console.log('Column by name:');
console.log(ageColumn);
ageColumn = d.columns.byIndex(3); 
console.log('Column by index:');
console.log(ageColumn);