// Different ways of iterating over column values have drastic impacts on performance

// Naive: 85ms
// Values iterator: 16ms
// Column instance: 14ms
// Column instance and length: 9ms
// Raw number array: 1ms

let t = grok.data.demo.demog(100000);
let sum = 0;

DG.time('Naive', () => {
  for (let i = 0; i < t.rowCount; i++)
    sum += t.columns.byName('age').get(i);
});

DG.time('Values iterator', () => {  
  let column = t.columns.byName('age');
  for (let v of column.values())
    sum += v;
});

DG.time('Column instance', () => {
  let column = t.columns.byName('age');
  for (let i = 0; i < t.rowCount; i++)
    sum += column.get(i);
});

DG.time('Column instance and length', () => {
  let column = t.columns.byName('age');
  let rowCount = t.rowCount;
  for (let i = 0; i < rowCount; i++)
    sum += column.get(i);
});

DG.time('Raw number array', () => {
  let array = t.columns.byName('age').getRawData();
  let rowCount = t.rowCount;
  for (let i = 0; i < rowCount; i++)
    sum += array[i];
});