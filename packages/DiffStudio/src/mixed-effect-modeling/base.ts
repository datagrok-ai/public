
// Seeded RNG (Mulberry32) for reproducibility
function mulberry32(seed: number) {
  let t = seed >>> 0;
  return function() {
    t += 0x6D2B79F5;
    let r = Math.imul(t ^ (t >>> 15), 1 | t);
    r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
    return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
  };
}

// Normal random using Box-Muller with seeded RNG
function normalRandom(rng: () => number, mean = 0, sd = 1) {
  // Box-Muller
  let u = 0; let v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  const z = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  return z * sd + mean;
}

// Helper: deep copy array of numbers
function cloneArray(a: number[]) {
  return a.slice();
}

// Pretty print table to console
function printResultsTable(table: { Parameter: string; TrueValue: number; Estimate: number; RelativeError: number }[]) {
  console.log('\nCOMPARISON OF TRUE AND ESTIMATED PARAMETERS:');
  console.table(table.map((r) => ({
    'Parameter': r.Parameter,
    'True Value': r.TrueValue,
    'Estimate': Number(r.Estimate.toFixed(6)),
    'Relative Error (%)': Number(r.RelativeError.toFixed(2)),
  })));
}

// -----------------------------
// STEP 1: DATA GENERATION
// -----------------------------
console.log('='.repeat(70));
console.log('STEP 1: DATA GENERATION');
console.log('='.repeat(70));

interface Observation {
  ID: number;
  Time: number;
  Concentration: number;
  CL_true: number;
  V_true: number;
  ka_true: number;
}

function generatePkData(options?: { nSubjects?: number; nTimepoints?: number; seed?: number }) {
  const nSubjects = options?.nSubjects ?? 10;
  const nTimepoints = options?.nTimepoints ?? 8;
  const seed = options?.seed ?? 42;

  const rng = mulberry32(seed);

  // Population parameters (fixed effects)
  const theta_CL = 5.0; // L/h
  const theta_V = 20.0; // L
  const theta_ka = 1.0; // 1/h

  // Inter-individual variability (standard deviations on log scale)
  const omega_CL = 0.3; // 30% on log-scale
  const omega_V = 0.2;
  const omega_ka = 0.4;

  // Residual error (proportional) — expressed as sd of log-error
  const sigma = 0.1; // log-scale sd

  const timePoints = [0.5, 1, 2, 4, 6, 8, 12, 24];
  const dose = 100;

  const data: Observation[] = [];

  for (let subjectId = 1; subjectId <= nSubjects; subjectId++) {
    // random effects (on log scale)
    const eta_CL = normalRandom(rng, 0, omega_CL);
    const eta_V = normalRandom(rng, 0, omega_V);
    const eta_ka = normalRandom(rng, 0, omega_ka);

    const CL_i = theta_CL * Math.exp(eta_CL);
    const V_i = theta_V * Math.exp(eta_V);
    const ka_i = theta_ka * Math.exp(eta_ka);

    for (const t of timePoints) {
      const ke = CL_i / V_i;
      let C_pred = 0;
      if (t > 0) {
        // handle case ka == ke carefully
        if (Math.abs(ka_i - ke) < 1e-12) {
          // limit when ka -> ke
          C_pred = (dose / V_i) * ka_i * t * Math.exp(-ke * t);
        } else
          C_pred = (dose / V_i) * (ka_i / (ka_i - ke)) * (Math.exp(-ke * t) - Math.exp(-ka_i * t));
      }

      const epsilon = normalRandom(rng, 0, sigma);
      const C_obs = Math.max(C_pred * Math.exp(epsilon), 1e-8);

      data.push({ID: subjectId, Time: t, Concentration: C_obs, CL_true: CL_i, V_true: V_i, ka_true: ka_i});
    }
  }

  console.log(`\n✓ Data generated for ${nSubjects} subjects`);
  console.log(`✓ Total observations: ${data.length}`);
  console.log('\nPopulation parameters (true values):');
  console.log(`  CL = ${theta_CL} L/h`);
  console.log(`  V  = ${theta_V} L`);
  console.log(`  ka = ${theta_ka} 1/h`);

  return {data, trueParams: {CL: theta_CL, V: theta_V, ka: theta_ka}};
}

const {data: dfData, trueParams} = generatePkData({nSubjects: 10, seed: 42});

// -----------------------------
// STEP 2: (Optional) Data inspection
// -----------------------------
console.log('\n'.repeat(1));
console.log('STEP 2: DATA QUICK INSPECTION');
console.log('Unique subjects:', Array.from(new Set(dfData.map((d) => d.ID))).length);
console.log('Unique time points:', Array.from(new Set(dfData.map((d) => d.Time))).length);

// -----------------------------
// STEP 3: MIXED-EFFECTS MODEL
// -----------------------------
console.log('\n' + '='.repeat(70));
console.log('STEP 3: BUILD MIXED-EFFECTS MODEL');
console.log('='.repeat(70));

class PKMixedEffectsModel {
  dose: number;

  constructor(dose = 100) {
    this.dose = dose;
  }

  predictConcentration(time: number, CL: number, V: number, ka: number) {
    const ke = CL / V;
    if (time <= 0 || ka <= ke)
      return 0.001; // stability floor


    if (Math.abs(ka - ke) < 1e-12) {
      const C = (this.dose / V) * ka * time * Math.exp(-ke * time);
      return Math.max(C, 0.001);
    }

    const C = (this.dose / V) * (ka / (ka - ke)) * (Math.exp(-ke * time) - Math.exp(-ka * time));
    return Math.max(C, 0.001);
  }

  // Negative log-likelihood (simplified)
  // params: [log(CL), log(V), log(ka), omega_CL, omega_V, omega_ka, log(sigma)]
  negativeLogLikelihood(params: number[], data: Observation[]) {
    const log_CL = params[0];
    const log_V = params[1];
    const log_ka = params[2];
    const theta_CL = Math.exp(log_CL);
    const theta_V = Math.exp(log_V);
    const theta_ka = Math.exp(log_ka);

    const omega_CL = Math.abs(params[3]);
    const omega_V = Math.abs(params[4]);
    const omega_ka = Math.abs(params[5]);

    const log_sigma = params[6];
    const sigma = Math.exp(log_sigma);

    let nll = 0;

    const uniqueSubjects = Array.from(new Set(data.map((d) => d.ID)));

    for (const subjectId of uniqueSubjects) {
      const subjectData = data.filter((d) => d.ID === subjectId);

      // Simplified EBE: use population values (no individual eta estimation here)
      const CL_i = theta_CL;
      const V_i = theta_V;
      const ka_i = theta_ka;

      for (const row of subjectData) {
        const C_pred = this.predictConcentration(row.Time, CL_i, V_i, ka_i);
        const C_obs = row.Concentration;
        if (C_obs > 0 && C_pred > 0) {
          const residual = Math.log(C_obs) - Math.log(C_pred);
          nll += 0.5 * Math.pow(residual / sigma, 2) + Math.log(sigma);
        } else {
          // penalize impossible values
          nll += 1e6;
        }
      }

      // penalty from random effects (we do not have subject-level etas here; just penalize omega sizes)
      nll += 0.5 * (omega_CL * omega_CL + omega_V * omega_V + omega_ka * omega_ka);
    }

    return nll;
  }
}

// -----------------------------
// Simple Nelder–Mead optimizer implementation (for demo purposes)
// This is a basic version and not production-grade.
// -----------------------------

interface NMOptions {
  maxIter?: number;
  tol?: number;
  alpha?: number; // reflection
  gamma?: number; // expansion
  rho?: number; // contraction
  sigma?: number; // shrink
}

function nelderMead(
  fn: (x: number[]) => number,
  x0: number[],
  options?: NMOptions,
): { x: number[]; fx: number; nit: number; success: boolean } {
  const n = x0.length;
  const maxIter = options?.maxIter ?? 1000;
  const tol = options?.tol ?? 1e-6;
  const alpha = options?.alpha ?? 1;
  const gamma = options?.gamma ?? 2;
  const rho = options?.rho ?? 0.5;
  const sigma = options?.sigma ?? 0.5;

  // initial simplex: x0 and x0 + e_i * scale
  const scale = 0.05;
  const simplex: number[][] = [cloneArray(x0)];
  for (let i = 0; i < n; i++) {
    const xi = cloneArray(x0);
    xi[i] = xi[i] + scale * (xi[i] === 0 ? 0.1 : xi[i]);
    simplex.push(xi);
  }

  let fs = simplex.map((x) => fn(x));

  let iter = 0;
  while (iter < maxIter) {
    // sort simplex by function values
    const idx = fs.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v).map((o) => o.i);
    const s = idx.map((i) => simplex[i]);
    const f = idx.map((i) => fs[i]);

    // check convergence
    const fMean = f.reduce((a, b) => a + b, 0) / f.length;
    const sq = f.map((fi) => Math.pow(fi - fMean, 2));
    const sd = Math.sqrt(sq.reduce((a, b) => a + b, 0) / f.length);
    if (sd < tol)
      return {x: s[0], fx: f[0], nit: iter, success: true};


    // centroid of all but worst
    const xBar = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++)
        xBar[i] += s[j][i];

      xBar[i] /= n;
    }

    // reflection
    const xR = xBar.map((xi, i) => xi + alpha * (xi - s[n][i]));
    const fR = fn(xR);

    if (fR < f[0]) {
      // expansion
      const xE = xBar.map((xi, i) => xi + gamma * (xR[i] - xi));
      const fE = fn(xE);
      if (fE < fR) {
        simplex[idx[n]] = xE;
        fs[idx[n]] = fE;
      } else {
        simplex[idx[n]] = xR;
        fs[idx[n]] = fR;
      }
    } else if (fR < f[n - 1]) {
      simplex[idx[n]] = xR;
      fs[idx[n]] = fR;
    } else {
      // contraction
      if (fR < f[n]) {
        const xC = xBar.map((xi, i) => xi + rho * (xR[i] - xi));
        const fC = fn(xC);
        if (fC <= fR) {
          simplex[idx[n]] = xC;
          fs[idx[n]] = fC;
        } else {
          // shrink
          for (let j = 1; j <= n; j++) {
            simplex[idx[j]] = s[0].map((x0i, i) => x0i + sigma * (s[j][i] - x0i));
            fs[idx[j]] = fn(simplex[idx[j]]);
          }
        }
      } else {
        const xC = xBar.map((xi, i) => xi + rho * (s[n][i] - xi));
        const fC = fn(xC);
        if (fC < f[n]) {
          simplex[idx[n]] = xC;
          fs[idx[n]] = fC;
        } else {
          // shrink
          for (let j = 1; j <= n; j++) {
            simplex[idx[j]] = s[0].map((x0i, i) => x0i + sigma * (s[j][i] - x0i));
            fs[idx[j]] = fn(simplex[idx[j]]);
          }
        }
      }
    }

    fs = simplex.map((x) => fn(x));
    iter += 1;
  }

  // return best found
  const bestIdx = fs.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v)[0].i;
  return {x: simplex[bestIdx], fx: fs[bestIdx], nit: iter, success: false};
}

// -----------------------------
// Fit model
// -----------------------------

console.log('\nInitializing model...');
const model = new PKMixedEffectsModel(100);

function fitModel(data: Observation[], verbose = true) {
  if (verbose) {
    console.log('Starting parameter estimation...');
    console.log('   Using (approximate) maximum likelihood with Nelder-Mead');
  }

  // Initial parameter guesses
  const initialParams = [Math.log(5.0), Math.log(20.0), Math.log(1.0), 0.3, 0.2, 0.4, Math.log(0.1)];

  // We'll not enforce strict bounds in this Nelder–Mead demo; instead we rely on penalties inside NLL.
  const objective = (x: number[]) => model.negativeLogLikelihood(x, data);

  const result = nelderMead(objective, initialParams, {maxIter: 2000, tol: 1e-5});

  if (!result.success)
    console.warn('\n Optimization did not converge to the requested tolerance.');


  const est = {
    CL: Math.exp(result.x[0]),
    V: Math.exp(result.x[1]),
    ka: Math.exp(result.x[2]),
    omega_CL: Math.abs(result.x[3]),
    omega_V: Math.abs(result.x[4]),
    omega_ka: Math.abs(result.x[5]),
    sigma: Math.exp(result.x[6]),
  };

  if (verbose) {
    console.log('\n✓ Estimation completed!');
    console.log(`   Iterations: ${result.nit}`);
    console.log(`   Final objective function value: ${result.fx.toFixed(6)}`);
  }

  return {estimates: est, optimizationResult: result};
}

const fitRes = fitModel(dfData, true);
const estimates = fitRes.estimates;

// -----------------------------
// STEP 4: RESULTS ANALYSIS
// -----------------------------
console.log('\n' + '='.repeat(70));
console.log('STEP 4: PARAMETER ESTIMATION RESULTS');
console.log('='.repeat(70));

const resultsTable = [
  {Parameter: 'CL (L/h)', TrueValue: trueParams.CL, Estimate: estimates.CL, RelativeError: Math.abs(estimates.CL - trueParams.CL) / trueParams.CL * 100},
  {Parameter: 'V (L)', TrueValue: trueParams.V, Estimate: estimates.V, RelativeError: Math.abs(estimates.V - trueParams.V) / trueParams.V * 100},
  {Parameter: 'ka (1/h)', TrueValue: trueParams.ka, Estimate: estimates.ka, RelativeError: Math.abs(estimates.ka - trueParams.ka) / trueParams.ka * 100},
  {Parameter: 'ω_CL', TrueValue: 0.3, Estimate: estimates.omega_CL, RelativeError: Math.abs(estimates.omega_CL - 0.3) / 0.3 * 100},
  {Parameter: 'ω_V', TrueValue: 0.2, Estimate: estimates.omega_V, RelativeError: Math.abs(estimates.omega_V - 0.2) / 0.2 * 100},
  {Parameter: 'ω_ka', TrueValue: 0.4, Estimate: estimates.omega_ka, RelativeError: Math.abs(estimates.omega_ka - 0.4) / 0.4 * 100},
  {Parameter: 'σ', TrueValue: 0.1, Estimate: estimates.sigma, RelativeError: Math.abs(estimates.sigma - 0.1) / 0.1 * 100},
];

printResultsTable(resultsTable);

// Compute predictions and residuals
const predictions: number[] = [];
for (const row of dfData) {
  const pred = model.predictConcentration(row.Time, estimates.CL, estimates.V, estimates.ka);
  predictions.push(pred);
}

// Attach predicted values to data (non-destructive)
const dfWithPred = dfData.map((r, i) => ({...r, Predicted: predictions[i]}));

// Residuals on log-scale
const residuals = dfWithPred.map((r) => Math.log(r.Concentration) - Math.log(r.Predicted));

// Basic diagnostic numbers
const maxObserved = Math.max(...dfWithPred.map((r) => r.Concentration));
const maxPredicted = Math.max(...dfWithPred.map((r) => r.Predicted));
const maxVal = Math.max(maxObserved, maxPredicted);

console.log('\nGoodness of fit:');
console.log(`  max observed: ${maxObserved.toFixed(6)}, max predicted: ${maxPredicted.toFixed(6)}, scale: ${maxVal.toFixed(6)}`);

console.log('\nResiduals summary (log-scale):');
const meanRes = residuals.reduce((a, b) => a + b, 0) / residuals.length;
const sdRes = Math.sqrt(residuals.map((r) => Math.pow(r - meanRes, 2)).reduce((a, b) => a + b, 0) / residuals.length);
console.log(`  mean: ${meanRes.toFixed(6)}, sd: ${sdRes.toFixed(6)}`);

// The arrays dfWithPred, predictions, residuals, and resultsTable are available for further analysis or plotting.

console.log('\nDone.');
