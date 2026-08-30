// Web Worker: brute-force grid search for Lotka-Volterra "max prey" optimization.
// Uses the MRT solver from diff-grok (same method as the main simulation).
import {mrt, ODEs} from 'diff-grok';

function maxPreyMRT(
  alpha: number, beta: number, delta: number, gamma: number,
  x0: number, y0: number, T: number,
): number {
  try {
    const odes: ODEs = {
      name: 'LV',
      arg: {name: 't', start: 0, finish: T, step: Math.max(0.1, T / 500)},
      initial: [x0, y0],
      func: (_t: number, y: Float64Array, out: Float64Array) => {
        out[0] = alpha * y[0] - beta * y[0] * y[1];
        out[1] = delta * y[0] * y[1] - gamma * y[1];
      },
      tolerance: 1e-4,
      solutionColNames: ['x', 'y'],
    };
    const sol = mrt(odes);
    const x = sol[1];
    let maxX = 0;
    for (let i = 0; i < x.length; i++) {
      if (!isFinite(x[i]) || x[i] > 1e8) return 0;
      if (x[i] > maxX) maxX = x[i];
    }
    return maxX;
  } catch {
    return 0;
  }
}

(self as any).onmessage = function(e: MessageEvent): void {
  const {alphaRange, betaRange, deltaRange, gammaRange, x0, y0, T} = e.data as {
    alphaRange: [number, number];
    betaRange: [number, number];
    deltaRange: [number, number];
    gammaRange: [number, number];
    x0: number; y0: number; T: number;
  };

  // 10% step = 10 intervals → 11 sample values per parameter
  const STEPS = 10;
  const aStep = (alphaRange[1] - alphaRange[0]) / STEPS;
  const bStep = (betaRange[1] - betaRange[0]) / STEPS;
  const dStep = (deltaRange[1] - deltaRange[0]) / STEPS;
  const gStep = (gammaRange[1] - gammaRange[0]) / STEPS;

  let bestAlpha = alphaRange[0];
  let bestBeta = betaRange[0];
  let bestDelta = deltaRange[0];
  let bestGamma = gammaRange[0];
  let bestMax = 0;

  const total = (STEPS + 1) ** 4;
  let done = 0;

  for (let ai = 0; ai <= STEPS; ai++) {
    const alpha = alphaRange[0] + ai * aStep;
    for (let bi = 0; bi <= STEPS; bi++) {
      const beta = betaRange[0] + bi * bStep;
      for (let di = 0; di <= STEPS; di++) {
        const delta = deltaRange[0] + di * dStep;
        for (let gi = 0; gi <= STEPS; gi++) {
          const gamma = gammaRange[0] + gi * gStep;
          const maxPrey = maxPreyMRT(alpha, beta, delta, gamma, x0, y0, T);
          if (maxPrey > bestMax) {
            bestMax = maxPrey;
            bestAlpha = alpha;
            bestBeta = beta;
            bestDelta = delta;
            bestGamma = gamma;
          }
          done++;
        }
      }
      // ~121 progress posts total (one per alpha×beta pair)
      (self as any).postMessage({type: 'progress', progress: done / total});
    }
  }

  (self as any).postMessage({
    type: 'result',
    alpha: bestAlpha,
    beta: bestBeta,
    delta: bestDelta,
    gamma: bestGamma,
    maxPrey: bestMax,
  });
};
