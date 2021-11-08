//name: Lotka-Volterra Model
//tags: model, NonlinearDynamics
//description: Model describes the dynamics of biological systems in which two species interact, one as a predator and the other as prey.
//language: javascript
//input: double x0 = 0.5 [Initial conditions for prey population]
//input: double y0 = 2 [Initial conditions for predator population]
//input: double alpha = 0.9 [Natural growth rate of preys in the absence of predation]
//input: double beta = 0.8 [Death rate per encounter of preys due to predation]
//input: double gamma = 0.9 [Natural death rate of predators in the absence of food (preys)]
//input: double sigma = 0.7 [Efficiency of turning eaten preys into predators]
//input: double h = 0.01 [Step height]
//input: double timeStart = 0
//input: double timeStop = 18
//output: dataframe df {viewer: Line Chart(x: "Time", multiAxis: "true", title: "Lotka-Volterra Model")}

function sumArrays(...arrays) {
  return Array.from({ length: 2 })
    .map((_, i) => arrays
      .map(xs => xs[i])
      .reduce((sum, x) => sum + x, 0)
    );
}

function f(r) {
  const fxd = r[0] * (alpha - beta * r[1]);
  const fyd = -r[1] * (gamma - sigma * r[0]);
  return new Float32Array([fxd, fyd]);
}

function rungeKuttaFourthOrderMethod(r, t, h) {
  const k1 = f(r).map((e) => e * h);
  const k2 = f(sumArrays(r, k1.map((e) => 0.5 * e))).map((e) => 2 * e * h);
  const k3 = f(sumArrays(r, k2.map((e) => 0.5 * e))).map((e) => 2 * e * h);
  const k4 = f(sumArrays(r, k3)).map((e) => e * h);
  return sumArrays(k1, k2, k3, k4);
}

let arange = function(start, stop, step) {
  return new Float32Array(Math.ceil((stop - start) / step))
    .map((_, i) => i * step);
};

const timePoints = arange(timeStart, timeStop, h);
let prey = new Float32Array(timePoints.length);
let predator = new Float32Array(timePoints.length);
let r = new Float32Array([x0, y0]);

timePoints.forEach((timePoint, index) => {
  prey[index] = r[0];
  predator[index] = r[1];
  r = sumArrays(r, rungeKuttaFourthOrderMethod(r, timePoint, h));
});

df = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('Time', timePoints),
  DG.Column.fromFloat32Array('Number of preys', prey),
  DG.Column.fromFloat32Array('Number of predators', predator)
]);