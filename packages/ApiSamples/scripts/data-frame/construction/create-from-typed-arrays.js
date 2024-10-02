//tags: DataFrame, Column, construction
// This is the fastest way of creating numerical columns
// The following example creates a table with 100 million rows and three columns
// See also: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays

let len = 100000000;
let ints = new Int32Array(len);
let floats = new Float32Array(len);
let floats64 = new Float64Array(len);

for (let i = 0; i < len; i++)
  ints[i] = i;

for (let i = 0; i < len; i++)
  floats[i] = i / 10;

for (let i = 0; i < len; i++)
  floats64[i] = i + i / len;

let table = DG.DataFrame.fromColumns([
  DG.Column.fromInt32Array('ints', ints),
  DG.Column.fromFloat32Array('floats', floats),
  DG.Column.fromFloat64Array('floats64', floats64)
]);

grok.shell.addTableView(table);
