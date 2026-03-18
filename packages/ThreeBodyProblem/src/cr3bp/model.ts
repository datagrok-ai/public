// CR3BP — ODE specification and mathematical utilities

import {ODEs} from 'diff-grok';

/** Parameters for the CR3BP simulation */
export interface CR3BPParams {
  mu: number;
  x0: number;
  y0: number;
  vx0: number;
  vy0: number;
  T: number;
}

/** Lagrange point data */
export interface LagrangePoint {
  name: string;
  x: number;
  y: number;
}

/** Zero-velocity curve grid data */
export interface ZVCGrid {
  xArr: Float64Array;
  yArr: Float64Array;
  U: Float64Array;
  cj0: number;
}

/** Creates the ODEs specification for the CR3BP */
export function createCR3BPODE(params: CR3BPParams): ODEs {
  const {mu, x0, y0, vx0, vy0, T} = params;

  return {
    name: 'CR3BP',
    arg: {name: 't', start: 0, finish: T, step: 0.001},
    initial: [x0, y0, vx0, vy0],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      const r1 = Math.sqrt((y[0] + mu) ** 2 + y[1] ** 2);
      const r2 = Math.sqrt((y[0] - 1 + mu) ** 2 + y[1] ** 2);
      const r1_3 = r1 ** 3;
      const r2_3 = r2 ** 3;
      out[0] = y[2]; // dx/dt = vx
      out[1] = y[3]; // dy/dt = vy
      out[2] = 2 * y[3] + y[0] - (1 - mu) * (y[0] + mu) / r1_3 - mu * (y[0] - 1 + mu) / r2_3;
      out[3] = -2 * y[2] + y[1] - (1 - mu) * y[1] / r1_3 - mu * y[1] / r2_3;
    },
    tolerance: 1e-8,
    solutionColNames: ['x', 'y', 'vx', 'vy'],
  };
}

/** Evaluates the CR3BP ODE right-hand side at a given state */
export function evalRHS(
  mu: number, x: number, y: number, vx: number, vy: number,
): [number, number, number, number] {
  const state = new Float64Array([x, y, vx, vy]);
  const out = new Float64Array(4);
  const r1 = Math.sqrt((state[0] + mu) ** 2 + state[1] ** 2);
  const r2 = Math.sqrt((state[0] - 1 + mu) ** 2 + state[1] ** 2);
  const r1_3 = r1 ** 3;
  const r2_3 = r2 ** 3;
  out[0] = state[2];
  out[1] = state[3];
  out[2] = 2 * state[3] + state[0] - (1 - mu) * (state[0] + mu) / r1_3 - mu * (state[0] - 1 + mu) / r2_3;
  out[3] = -2 * state[2] + state[1] - (1 - mu) * state[1] / r1_3 - mu * state[1] / r2_3;
  return [out[0], out[1], out[2], out[3]];
}

/** Computes the effective potential U(x, y) */
export function effectivePotential(x: number, y: number, mu: number): number {
  const r1 = Math.sqrt((x + mu) ** 2 + y ** 2);
  const r2 = Math.sqrt((x - 1 + mu) ** 2 + y ** 2);
  return 0.5 * (x ** 2 + y ** 2) + (1 - mu) / r1 + mu / r2;
}

/** Computes the Jacobi constant */
export function jacobiConstant(x: number, y: number, vx: number, vy: number, mu: number): number {
  const U = effectivePotential(x, y, mu);
  return 2 * U - (vx ** 2 + vy ** 2);
}

/** Finds a collinear Lagrange point via Newton's method */
function findCollinearLagrange(mu: number, x0: number, maxIter: number = 100, tol: number = 1e-12): number {
  let x = x0;

  for (let i = 0; i < maxIter; i++) {
    const r1 = Math.abs(x + mu);
    const r2 = Math.abs(x - 1 + mu);
    if (r1 < 1e-15 || r2 < 1e-15) break;

    const r1_3 = r1 ** 3;
    const r2_3 = r2 ** 3;
    const f = x - (1 - mu) * (x + mu) / r1_3 - mu * (x - 1 + mu) / r2_3;

    // Numerical derivative via central difference
    const h = 1e-8;
    const xp = x + h;
    const xm = x - h;
    const r1p = Math.abs(xp + mu);
    const r2p = Math.abs(xp - 1 + mu);
    const r1m = Math.abs(xm + mu);
    const r2m = Math.abs(xm - 1 + mu);
    const fp = xp - (1 - mu) * (xp + mu) / (r1p ** 3) - mu * (xp - 1 + mu) / (r2p ** 3);
    const fm = xm - (1 - mu) * (xm + mu) / (r1m ** 3) - mu * (xm - 1 + mu) / (r2m ** 3);
    const df = (fp - fm) / (2 * h);

    if (Math.abs(df) < 1e-15) break;

    const dx = -f / df;
    x += dx;
    if (Math.abs(dx) < tol) break;
  }

  return x;
}

/** Computes all 5 Lagrange points for given mass parameter */
export function computeLagrangePoints(mu: number): LagrangePoint[] {
  // Hill's approximation for initial guesses
  const cbrtMu3 = Math.cbrt(mu / 3);

  // L1: between bodies
  const xL1 = findCollinearLagrange(mu, 1 - mu - cbrtMu3);

  // L2: beyond smaller body
  const xL2 = findCollinearLagrange(mu, 1 - mu + cbrtMu3);

  // L3: beyond larger body
  const xL3 = findCollinearLagrange(mu, -1 - 5 * mu / 12);

  // L4, L5: triangular points (analytical)
  const xL4 = 0.5 - mu;
  const yL4 = Math.sqrt(3) / 2;

  return [
    {name: 'L1', x: xL1, y: 0},
    {name: 'L2', x: xL2, y: 0},
    {name: 'L3', x: xL3, y: 0},
    {name: 'L4', x: xL4, y: yL4},
    {name: 'L5', x: xL4, y: -yL4},
  ];
}

/** Computes the zero-velocity curve grid */
export function computeZVCGrid(mu: number, cj0: number, gridSize: number = 200): ZVCGrid {
  const xMin = -1.5;
  const xMax = 1.5;
  const yMin = -1.5;
  const yMax = 1.5;

  // First pass: count forbidden points
  let count = 0;
  for (let i = 0; i < gridSize; i++) {
    for (let j = 0; j < gridSize; j++) {
      const gx = xMin + (xMax - xMin) * i / (gridSize - 1);
      const gy = yMin + (yMax - yMin) * j / (gridSize - 1);
      const U2 = 2 * effectivePotential(gx, gy, mu);
      if (U2 > cj0) count++;
    }
  }

  // Second pass: fill arrays
  const xArr = new Float64Array(count);
  const yArr = new Float64Array(count);
  const UArr = new Float64Array(count);
  let idx = 0;

  for (let i = 0; i < gridSize; i++) {
    for (let j = 0; j < gridSize; j++) {
      const gx = xMin + (xMax - xMin) * i / (gridSize - 1);
      const gy = yMin + (yMax - yMin) * j / (gridSize - 1);
      const U2 = 2 * effectivePotential(gx, gy, mu);
      if (U2 > cj0) {
        xArr[idx] = gx;
        yArr[idx] = gy;
        UArr[idx] = U2;
        idx++;
      }
    }
  }

  return {xArr, yArr, U: UArr, cj0};
}
