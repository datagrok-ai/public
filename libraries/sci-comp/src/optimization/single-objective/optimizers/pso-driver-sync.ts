/* eslint-disable */
// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./pso-driver.ts
import {OptimizationResult} from '../types';
import {PSOSettings} from './pso';
import {PSODeps, Particle, RNG, defaultLower, defaultUpper, mulberry32} from './pso-driver';

function initSwarm(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  size: number,
  lo: Float64Array,
  hi: Float64Array,
  vMax: Float64Array,
  rng: RNG,
): Particle[] {
  const n = x0.length;
  const swarm: Particle[] = [];

  // First particle starts at x0
  const p0 = new Float64Array(n);
  p0.set(x0);
  const val0 = fn(p0);
  const bp0 = new Float64Array(n);
  bp0.set(p0);
  swarm.push({
    position: p0,
    velocity: new Float64Array(n),
    value: val0,
    bestPosition: bp0,
    bestValue: val0,
  });

  // Rest — random within search range
  for (let k = 1; k < size; k++) {
    const pos = new Float64Array(n);
    const vel = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      pos[i] = lo[i] + rng() * (hi[i] - lo[i]);
      vel[i] = -vMax[i] + rng() * 2 * vMax[i];
    }

    const val = fn(pos);
    const bestPos = new Float64Array(n);
    bestPos.set(pos);

    swarm.push({
      position: pos,
      velocity: vel,
      value: val,
      bestPosition: bestPos,
      bestValue: val,
    });
  }

  return swarm;
}

export function runPSOSync(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  s: PSOSettings,
  deps: PSODeps,
): OptimizationResult {
  const n = x0.length;
  const swarmSize = s.swarmSize!;
  const maxIter = s.maxIterations!;
  const tol = s.tolerance!;
  const w = s.inertia!;
  const c1 = s.cognitive!;
  const c2 = s.social!;
  const rng: RNG = s.seed != null ? mulberry32(s.seed) : Math.random;

  const rangeLo = s.searchRange?.lower ?? defaultLower(x0);
  const rangeHi = s.searchRange?.upper ?? defaultUpper(x0);
  const vMax = new Float64Array(n);
  for (let i = 0; i < n; i++)
    vMax[i] = s.maxVelocityFraction! * (rangeHi[i] - rangeLo[i]);

  // --- Initialize swarm ---
  const swarm = initSwarm(fn, x0, swarmSize, rangeLo, rangeHi, vMax, rng);

  // Global best
  const gBest = new Float64Array(n);
  gBest.set(swarm[0].bestPosition);
  let gBestVal = swarm[0].bestValue;
  for (const p of swarm) {
    if (p.bestValue < gBestVal) {
      gBestVal = p.bestValue;
      gBest.set(p.bestPosition);
    }
  }

  const costHistory = new Float64Array(maxIter);
  let costLen = 0;
  let iteration = 0;
  let prevBest = Infinity;
  let noImprovement = 0;
  let converged = false;

  while (iteration < maxIter) {
    costHistory[costLen++] = gBestVal;

    // --- Convergence ---
    if (iteration > 0 && prevBest - gBestVal > tol)
      noImprovement = 0;
    else if (iteration > 0) {
      noImprovement++;
      if (noImprovement >= s.noImprovementLimit!) {converged = true; break;}
    }
    prevBest = gBestVal;

    // --- Callback ---
    if (deps.notify(s.onIteration, {
      iteration,
      bestValue: gBestVal,
      bestPoint: gBest,
    }))
      break;

    // --- Update particles ---
    for (const p of swarm) {
      for (let i = 0; i < n; i++) {
        const r1 = rng();
        const r2 = rng();

        p.velocity[i] =
          w * p.velocity[i] +
          c1 * r1 * (p.bestPosition[i] - p.position[i]) +
          c2 * r2 * (gBest[i] - p.position[i]);

        if (p.velocity[i] > vMax[i]) p.velocity[i] = vMax[i];
        else if (p.velocity[i] < -vMax[i]) p.velocity[i] = -vMax[i];

        p.position[i] += p.velocity[i];
      }

      // Evaluate
      p.value = fn(p.position);

      // Update personal best
      if (p.value < p.bestValue) {
        p.bestValue = p.value;
        p.bestPosition.set(p.position);

        // Update global best
        if (p.value < gBestVal) {
          gBestVal = p.value;
          gBest.set(p.position);
        }
      }
    }

    iteration++;
  }

  return {
    point: gBest,
    value: gBestVal,
    iterations: iteration,
    converged,
    costHistory: costHistory.subarray(0, costLen),
  };
}
