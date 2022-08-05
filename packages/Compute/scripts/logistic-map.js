//name: Logistic Map
//tags: simulation
//description: Archetypal example of how complex, chaotic behaviour can arise from very simple non-linear dynamical equations.
//language: javascript
//input: int iterations = 1000 [Number of iterations per point]
//input: double x0 = 0.5 [Initial conditions]
//input: double spacing = 0.001 [Spacing between points on domain (r-axis)]
//input: double res = 8 [Largest n-cycle visible]
//output: dataframe df {viewer: Scatter plot(title: "Logistic map")}
//meta.domain: Nonlinear dynamics

function getRandomNumberInRange(min, max) {
  return Math.random() * (max - min) + min;
}

function iterate(n, x, r) {
  for (let i = 1; i < n; i++)
    x = x * r * (1 - x);
  return x;
}

const start = Math.floor(1 / spacing);
const end = Math.ceil(4 / spacing);
const leftBoundary = iterations - res / 2;
const rightBoundary = iterations + res / 2;
let rs = new Float32Array(end - start);
let xs = new Float32Array(end - start);

for (let i = start; i < end; i++) {
  rs[i - start] = i * spacing;
  xs[i - start] = iterate(
    getRandomNumberInRange(leftBoundary, rightBoundary),
    x0,
    rs[i - start]
  );
}

df = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('R', rs),
  DG.Column.fromFloat32Array('X', xs)
]);