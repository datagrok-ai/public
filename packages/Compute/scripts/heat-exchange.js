//name: Heat Exchange
//language: javascript
//tags: simulation
//input: int len = 1000
//input: double k = 1.2
//output: dataframe simulation {viewer: Line Chart(multiAxis: "true")}
//output: double finalTemperature
//output: double finalSaturation
//output: double finalConcentration
//meta.domain: Manufacturing
//meta.modality: Small molecule

function sigmoid(z) {
  return 1 / (1 + Math.exp(-z));
}

let xs = new Float32Array(len);
let conc = new Float32Array(len);
let temp = new Float32Array(len);
let sat = new Float32Array(len);

for (let i = 0; i < len; i++) {
  let x = -5 + (i / len) * 10;
  xs[i] = x;
  temp[i] = 70 + Math.sin(x);
  sat[i] = 10 - Math.floor(10 * i / len);
  conc[i] = (temp[i] / 70)  + sigmoid(x);
}

finalTemperature = temp[len -1];
finalSaturation = sat[len - 1];
finalConcentration = conc[len - 1];
simulation = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('time', xs),
  DG.Column.fromFloat32Array('concentration', conc),
  DG.Column.fromFloat32Array('temp', temp),
  DG.Column.fromFloat32Array('saturation', sat),
]);