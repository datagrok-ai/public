// Problems for testing solvers' performance

/** Robertson chemical reaction, updated version (see https://archimede.uniba.it/~testset/problems/rober.php) */
const robertson = {
  name: 'Robertson',
  arg: {name: 't', start: 0, finish: 10e11, step: 10e9},
  initial: [1, 0, 0],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    output[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    output[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1]**2;
    output[2] = 3e7 * y[1]**2;
  },
  tolerance: 1e-7,
  solutionColNames: ['A', 'B', 'C'],
};

/** High Irradiance Responses of photomorphogenesis (see https://archimede.uniba.it/~testset/problems/hires.php) */
const hires = {
  name: 'HIRES',
  arg: {name: 't', start: 0, finish: 321.8122, step: 0.01},
  initial: [1, 0, 0, 0, 0, 0, 0, 0.0057],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    // extract function values
    const y1 = y[0];
    const y2 = y[1];
    const y3 = y[2];
    const y4 = y[3];
    const y5 = y[4];
    const y6 = y[5];
    const y7 = y[6];
    const y8 = y[7];

    // compute output
    output[0] = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007;
    output[1] = 1.71 * y1 - 8.75 * y2;
    output[2] = -10.03 * y3 + 0.43 * y4 + 0.035 * y5;
    output[3] = 8.32 * y2 + 1.71 * y3 - 1.12 * y4;
    output[4] = -1.745 * y5 + 0.43 * y6 + 0.43 * y7;
    output[5] = -280 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7;
    output[6] = 280 * y6 * y8 - 1.81 * y7;
    output[7] = -280 * y6 * y8 + 1.81 * y7;
  },
  tolerance: 1e-10,
  solutionColNames: ['y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8'],
};

/** Van der Pol oscillator (see https://archimede.uniba.it/~testset/problems/vdpol.php) */
const vanDerPol = {
  name: 'van der Pol',
  arg: {name: 't', start: 0, finish: 2000, step: 0.1},
  initial: [-1, 1],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    output[0] = y[1];
    output[1] = -y[0] + 1000 * (1 - y[0] * y[0]) * y[1];
  },
  tolerance: 1e-12,
  solutionColNames: ['x1', 'x2'],
};

export const problems = new Map([
  [robertson.name, robertson],
  [hires.name, hires],
  [vanDerPol.name, vanDerPol],
]);
