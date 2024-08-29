/* Problems for testing solvers' performance.
   Source: https://archimede.uniba.it/~testset/testsetivpsolvers/?page_id=26#ODE
*/

/** Robertson chemical reaction, updated version (see https://archimede.uniba.it/~testset/problems/rober.php) */
const robertson = {
  name: 'Robertson',
  arg: {name: 't', start: 0, finish: 10e11, step: 2.5e7},
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
}; // hires

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

/** The OREGO model (see https://archimede.uniba.it/~testset/report/orego.pdf) */
const orego = {
  name: 'OREGO',
  arg: {name: 't', start: 0, finish: 360, step: 0.01},
  initial: [1, 2, 3],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    // extract function values
    const y1 = y[0];
    const y2 = y[1];
    const y3 = y[2];

    // compute output
    output[0] = 77.27 * (y2 - y1 * y2 + y1 - 0.000008375 * y1 * y1 );
    output[1] = 1 / 77.27 * (-y2 - y1 * y2 + y3);
    output[2] = 0.161 * (y1 - y3);
  },
  tolerance: 1e-8,
  solutionColNames: ['y1', 'y2', 'y3'],
};

/** Kintetic constants for the E5 model */
enum E5 {
  K1 = 7.89e-10,
  K2 = 1.13e9,
  K3 = 1.1e7,
  K4 = 1.13e3,
};

/** The E5 model (chemical pyrolysis: https://archimede.uniba.it/~testset/report/e5.pdf) */
const e5 = {
  name: 'E5',
  arg: {name: 't', start: 0, finish: 1e13, step: 2.5e8},
  initial: [0.00176, 0, 0, 0],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    // extract function values
    const y1 = y[0];
    const y2 = y[1];
    const y3 = y[2];
    const y4 = y[3];

    // compute output
    output[0] = -E5.K1 * y1 - E5.K3 * y1 * y3;
    output[1] = E5.K1 * y1 - E5.K2 * y2 * y3;
    output[2] = E5.K1 * y1 - E5.K2 * y2 * y3 - E5.K3 * y1 * y3 + E5.K4 * y4;
    output[3] = E5.K3 * y1 * y3 - E5.K4 * y4;
  },
  tolerance: 1e-6,
  solutionColNames: ['y1', 'y2', 'y3', 'y4'],
};

/** Kintetic constants for the Pollution model */
enum POL {
  K1 = 0.35,
  K2 = 26.6,
  K3 = 1.23e4,
  K4 = 8.6e-4,
  K5 = 8.2e-4,
  K6 = 1.5e4,
  K7 = 1.3e-4,
  K8 = 2.4e4,
  K9 = 1.65e4,
  K10 = 9e3,
  K11 = 0.022,
  K12 = 1.2e4,
  K13 = 1.88,
  K14 = 1.63e4,
  K15 = 4.8e6,
  K16 = 3.5e-4,
  K17 = 0.0175,
  K18 = 1e8,
  K19 = 4.44e11,
  K20 = 1240,
  K21 = 2.1,
  K22 = 5.78,
  K23 = 0.0474,
  K24 = 1780,
  K25 = 3.12,
}; // POL

/** The chemical reaction part of the air pollution model (https://archimede.uniba.it/~testset/report/pollu.pdf) */
const pollution = {
  name: 'Pollution',
  arg: {name: 't', start: 0, finish: 60, step: 0.002},
  initial: [0, 0.2, 0, 0.04, 0, 0, 0.1, 0.3, 0.01, 0, 0, 0, 0, 0, 0, 0, 0.007, 0, 0, 0],
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
    const y9 = y[8];
    const y10 = y[9];
    const y11 = y[10];
    const y12 = y[11];
    const y13 = y[12];
    const y14 = y[13];
    const y15 = y[14];
    const y16 = y[15];
    const y17 = y[16];
    const y18 = y[17];
    const y19 = y[18];
    const y20 = y[19];

    // evaluate expressions
    const r1 = POL.K1 * y1;
    const r2 = POL.K2 * y2 * y4;
    const r3 = POL.K3 * y5 * y2;
    const r4 = POL.K4 * y7;
    const r5 = POL.K5 * y7;
    const r6 = POL.K6 * y7 * y6;
    const r7 = POL.K7 * y9;
    const r8 = POL.K8 * y9 * y6;
    const r9 = POL.K9 * y11 * y2;
    const r10 = POL.K10 * y11 * y1;
    const r11 = POL.K11 * y13;
    const r12 = POL.K12 * y10 * y2;
    const r13 = POL.K13 * y14;
    const r14 = POL.K14 * y1 * y6;
    const r15 = POL.K15 * y3;
    const r16 = POL.K16 * y4;
    const r17 = POL.K17 * y4;
    const r18 = POL.K18 * y16;
    const r19 = POL.K19 * y16;
    const r20 = POL.K20 * y17 * y6;
    const r21 = POL.K21 * y19;
    const r22 = POL.K22 * y19;
    const r23 = POL.K23 * y1 * y4;
    const r24 = POL.K24 * y19 * y1;
    const r25 = POL.K25 * y20;

    // compute output
    output[0] = -(r1 + r10 + r14 + r23 + r24) + (r2 + r3 + r9 + r11 + r12 + r22 + r25);
    output[1] = -r2 - r3 - r9 - r12 + r1 + r21;
    output[2] = -r15 + r1 + r17 + r19 + r22;
    output[3] = -r2 - r16 - r17 - r23 + r15;
    output[4] = -r3 +2 * r4 + r6 +r7 +r13 + r20;
    output[5] = -r6 - r8 - r14 - r20 + r3 + 2 * r18;
    output[6] = -r4 - r5 - r6 + r13;
    output[7] = r4 + r5 + r6 + r7;
    output[8] = -r7 - r8;
    output[9] = -r12 +r7 + r9;
    output[10] = -r9 - r10 + r8 + r11;
    output[11] = r9;
    output[12] = -r11 + r10;
    output[13] = -r13 + r12;
    output[14] = r14;
    output[15] = -r18 - r19 + r16;
    output[16] = -r20;
    output[17] = r20;
    output[18] = -r21 - r22 - r24 + r23 + r25;
    output[19] = -r25 + r24;
  },
  tolerance: 1e-6,
  solutionColNames: ['y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10', 'y11',
    'y12', 'y13', 'y14', 'y15', 'y16', 'y17', 'y18', 'y19', 'y20'],
}; // pollution

/** Problems for testing solvers' performance */
export const problems = new Map([
  [robertson.name, robertson],
  [hires.name, hires],
  [vanDerPol.name, vanDerPol],
  [orego.name, orego],
  [e5.name, e5],
  [pollution.name, pollution],
]);
