// This is the fastest way of creating numerical columns
// The following example creates a table with 100 million rows and two columns
// See also: https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays

let n = 100;
let len = n * n;
let x = new Float32Array(len);
let y = new Float32Array(len);
let z = new Float32Array(len);

for (let i = 0; i < len; i++) {
  x[i] = Math.floor(i / n) - n / 2;
  y[i] = i % n - n / 2;
  z[i] = x[i] * y[i];
}

let table = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('x', x),
  DG.Column.fromFloat32Array('y', y),
  DG.Column.fromFloat32Array('z', z),
]);

let tv = grok.shell.addTableView(table);
tv.addViewer('SurfacePlot');