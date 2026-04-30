//name: Lorenz Attractor
//meta.domain: Nonlinear dynamics
//description: Set of chaotic solutions of the system for weather prediction. It illustrates the phenomenon known as the Butterfly effect or (more technically) sensitive dependence on initial conditions
//language: javascript
//input: int iterations = 1000 {caption: Iterations number}
//input: double dt = 0.01 {caption: Step per iteration}
//input: double x0 = 0 [Initial conditions for rate of convective motion - i.e. how fast the rolls are rotating]
//input: double y0 = 1 [Initial conditions for temperature difference between the ascending and descending currents]
//input: double z0 = 1.05 [Initial conditions for distortion (from linearity) of the vertical temperature profile]
//output: dataframe df {caption: Lorenz attractor; viewer: 3d Scatter plot()}
//editor: Compute:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true

function lorenzEquations(x, y, z, s=10, r=28, b=2.667) {
  return {
    x: s * (y - x),
    y: r * x - y - x * z,
    z: x * y - b * z
  };
}

let xs = new Float32Array(iterations + 1);
let ys = new Float32Array(iterations + 1);
let zs = new Float32Array(iterations + 1);

xs[0] = x0;
ys[0] = y0;
zs[0] = z0;

for (let i = 0; i < iterations; i++) {
  let point = lorenzEquations(xs[i], ys[i], zs[i]);
  xs[i + 1] = xs[i] + point.x * dt;
  ys[i + 1] = ys[i] + point.y * dt;
  zs[i + 1] = zs[i] + point.z * dt;
}

df = DG.DataFrame.fromColumns([
  DG.Column.fromFloat32Array('X', xs),
  DG.Column.fromFloat32Array('Y', ys),
  DG.Column.fromFloat32Array('Z', zs)
]);
