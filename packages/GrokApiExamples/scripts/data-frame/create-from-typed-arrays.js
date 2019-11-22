// This is the fastest way of creating numerical columns
// The following example creates a table with 100 million rows and two columns
// See also: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays

let len = 100000000;
let ints = new Int32Array(len);
let floats = new Float32Array(len);

for (let i = 0; i < len; i++)
    ints[i] = i;

for (let i = 0; i < len; i++)
    floats[i] = i / 10;

let table = DataFrame.fromColumns([
   Column.fromInt32Array('ints', ints),
   Column.fromFloat32Array('floats', floats)
]);

gr.addTableView(table);
