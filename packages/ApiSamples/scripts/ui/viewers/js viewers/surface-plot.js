const n = 100;
const len = n * n;
const x = new Float32Array(len);
const y = new Float32Array(len);
const z = new Float32Array(len);

for (let i = 0; i < len; i++) {
  x[i] = Math.floor(i / n) - n / 2;
  y[i] = i % n - n / 2;
  z[i] = x[i] * y[i];
}

const table = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('x', x),
  DG.Column.fromFloat32Array('y', y),
  DG.Column.fromFloat32Array('z', z),
]);

let view = grok.shell.addTableView(table);
view.addViewer('surface plot');