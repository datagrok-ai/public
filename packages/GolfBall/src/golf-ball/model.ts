// Golf Ball Flight — ODE specification

import {ODEs} from 'diff-grok';

/** Parameters for the golf ball flight ODE system */
export interface GolfBallParams {
  v0: number;
  theta: number;
  m: number;
  d: number;
  Cd: number;
  rho: number;
  g: number;
}

/** Creates the ODEs specification for the golf ball flight model */
export function createGolfBallODE(params: GolfBallParams): ODEs {
  const {v0, theta, m, d, Cd, rho, g} = params;
  const thetaRad = theta * Math.PI / 180;
  const A = Math.PI * (d / 2) ** 2;
  const vx0 = v0 * Math.cos(thetaRad);
  const vy0 = v0 * Math.sin(thetaRad);

  return {
    name: 'GolfBallFlight',
    arg: {name: 't', start: 0, finish: 30, step: 0.01},
    initial: [0, 0, vx0, vy0],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      const vx = y[2];
      const vy = y[3];
      const speed = Math.sqrt(vx * vx + vy * vy);
      const drag = speed > 0 ? Cd * rho * A * speed / (2 * m) : 0;
      out[0] = vx;             // dx/dt
      out[1] = vy;             // dy/dt
      out[2] = -drag * vx;     // dvx/dt
      out[3] = -g - drag * vy; // dvy/dt
    },
    tolerance: 1e-6,
    solutionColNames: ['x(t)', 'y(t)', 'vx(t)', 'vy(t)'],
  };
}
