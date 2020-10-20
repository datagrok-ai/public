// Uses ColumnList as iterable and prints names of columns

let demog = grok.data.demo.demog();

for (let column of demog.columns)
    console.log(column.name);