let table = grok.data.testData('random walk', 100000, 2);

measure = function(f) {
  let t0 = performance.now();
  f();
  let t1 = performance.now();
  return t1 - t0;
};

let bitsetIteration = measure(function() {
  for (let i = 0; i < table.rowCount; i++)
    table.filter.get(i);
});

let bitsetIterationFast = measure(function() {
  let indexes = table.filter.getSelectedIndexes();
  let count = 0;
  for (let i = 0; i < indexes.length; i++)
    count++;
  console.log(`${count}`);
});

let valuesIteration = measure(function() {
  let col = table.columns.byIndex(0);
  for (let i = 0; i < table.rowCount; i++)
    col.get(i);
});

let valuesIterationFast = measure(function() {
  let sum = 0;
  let col = table.columns.byIndex(0);
  let data = col.getRawData();
  for (let i = 0; i < data.length; i++)
    sum += data[i];
});

ui.dialog()
  .add(`bitset: ${bitsetIteration}`)
  .add(`bitset fast: ${bitsetIterationFast}`)
  .add(`values: ${valuesIteration}`)
  .add(`values fast: ${valuesIterationFast}`)
  .show();
