/**
 * One global additive (Free-Wilson) fit over every role column of a decomposition at once,
 * `y ≈ μ + Σ_r θ_r(value)`, so one run puts every value of every role on a single scale.
 *
 * Legitimate only where a role value means the same thing in every row — a named fragment column.
 * A substituent label discovered by fragmentation is local to its own series, and an offset pooled
 * across series would then be averaging unrelated quantities.
 */

/** The ridge, applied as `acc / (count + prior)`: one notional compound at the library average added
 *  to every level. In units of count rather than activity, so it means the same in any assay units.
 *  Unpenalised, a design whose roles nearly alias each other grows offsets in the aliased direction
 *  that cancel in the prediction, and the per-role spread then reports noise as signal — a shuffled
 *  activity column produces larger spreads than the real one. */
const ROLE_FIT_PRIOR = 1;
/** Offsets are read at two decimals, so this is four orders below anything on screen. */
const ROLE_FIT_TOLERANCE = 1e-6;
/** A near-aliased design needs hundreds of sweeps. Stopping short leaves every role shrunk toward
 *  zero by a different amount, which reorders the roles without anything failing. */
export const ROLE_FIT_MAX_SWEEPS = 2000;
/** Floor under the size-aware cap: a design wide enough to hit it is one where even a converged fit
 *  would be too slow to run on the main thread, and a short run that reports itself unconverged
 *  suppresses the ranking rather than mis-stating it. */
const ROLE_FIT_MIN_SWEEPS = 50;
/** Sweeps scale with collinearity, not rows, so bounding rows alone refuses large easy tables and
 *  admits small expensive ones. */
const ROLE_FIT_WORK_BUDGET = 40_000_000;
const ROLE_CV_FOLDS = 5;
/** Below this many held-out points an R² is noise rather than a quality. */
const ROLE_MIN_CV_POINTS = 4;

export interface RoleDesign {
  /** Role column names, core first. */
  names: string[];
  /** values[r][k] — role r's value on observation k. Parallel arrays, one entry per measured cell. */
  values: string[][];
  activity: number[];
  /** Source-table row of each observation; the cross-validation folds key on it. */
  molIdx: number[];
  /** Observations a level needs before its offset is readable. */
  minSupport: number;
  /** Which end of the activity column is the better one. The fit itself is direction-free — it is
   *  in the column's own units — but `levels` is ordered best-first, and on a raw IC50 column the
   *  best level is the lowest. */
  higherIsBetter: boolean;
}

export interface RoleLevel {
  value: string;
  coef: number;
  n: number;
}

/** Two clusters of values separated by more than the fit can resolve. */
interface RoleSplit {
  gap: number;
  hiCount: number; hiMean: number; hiN: number;
  loCount: number; loMean: number; loN: number;
}

export interface RoleSummary {
  name: string;
  /** Levels carrying at least `minSupport` observations, best first. */
  levels: RoleLevel[];
  /** Levels carrying any observation at all, and how many of those are below the support floor. */
  fitted: number;
  thin: number;
  /** Count-weighted sd of this role's offsets with estimation noise subtracted. */
  spread: number;
  /** Correlation of this role's offsets between two halves of the table; null below three comparable
   *  levels, where a correlation is ±1 whatever the data and would read as evidence either way. */
  repeat: number | null;
  split: RoleSplit | null;
}

export interface RoleFit {
  roles: RoleSummary[];
  /** Count-weighted mean activity of the fitted compounds — what every offset is read against. */
  mean: number;
  compounds: number;
  /** Observations outside the largest connected component, which are not fitted. */
  dropped: number;
  sweeps: number;
  converged: boolean;
  residualSd: number;
  cvR2: number | null;
  cvRmse: number | null;
}

/** One solved backfit over the observations named by `rows`, keyed by the caller's level codes so a
 *  refit on a subset stays comparable with the full fit level for level. */
interface Backfit {
  mu: number;
  theta: Float64Array[];
  counts: Int32Array[];
  fitted: Float64Array;
  sweeps: number;
  converged: boolean;
}

/**
 * Averaging each role's margins in a single pass solves this model only when every combination is
 * measured the same number of times, which a library never is.
 */
function backfit(codes: Int32Array[], levelCounts: number[], y: Float64Array,
  rows: Int32Array): Backfit {
  const roleCount = codes.length;
  const m = rows.length;
  const counts = levelCounts.map((levels) => new Int32Array(levels));
  let sum = 0;
  for (let i = 0; i < m; i++) {
    const k = rows[i];
    sum += y[k];
    for (let r = 0; r < roleCount; r++)
      counts[r][codes[r][k]]++;
  }

  let mu = m ? sum / m : 0;
  const theta = levelCounts.map((levels) => new Float64Array(levels));
  const next = levelCounts.map((levels) => new Float64Array(levels));
  const fitted = new Float64Array(y.length).fill(mu);
  let sweeps = 0;
  let converged = false;
  const maxSweeps = Math.max(ROLE_FIT_MIN_SWEEPS,
    Math.min(ROLE_FIT_MAX_SWEEPS, Math.floor(ROLE_FIT_WORK_BUDGET / Math.max(1, m * roleCount))));
  while (!converged && sweeps < maxSweeps) {
    sweeps++;
    let moved = 0;
    for (let r = 0; r < roleCount; r++) {
      const code = codes[r];
      const n = counts[r];
      const th = theta[r];
      const acc = next[r];
      acc.fill(0);
      for (let i = 0; i < m; i++) {
        const k = rows[i];
        acc[code[k]] += y[k] - fitted[k] + th[code[k]];
      }
      let weighted = 0;
      let total = 0;
      for (let v = 0; v < acc.length; v++) {
        if (n[v] === 0)
          continue;
        acc[v] /= n[v] + ROLE_FIT_PRIOR;
        weighted += n[v] * acc[v];
        total += n[v];
      }
      // The prior pulls the whole role off centre by a constant, so the shift is taken out and folded
      // into the intercept: without it "against the library average" is simply false. It is also why
      // the movement below has to be measured after the re-centring — measured before, it plateaus at
      // that constant and the fit never reports convergence.
      const shift = total ? weighted / total : 0;
      mu += shift;
      for (let v = 0; v < acc.length; v++) {
        if (n[v] === 0)
          continue;
        acc[v] -= shift;
        moved = Math.max(moved, Math.abs(acc[v] - th[v]));
      }
      for (let i = 0; i < m; i++) {
        const k = rows[i];
        fitted[k] += shift + acc[code[k]] - th[code[k]];
      }
      th.set(acc);
    }
    // Every role in the sweep, not just the last one updated: with more than two roles the last can
    // be still while an earlier one is moving.
    converged = moved < ROLE_FIT_TOLERANCE;
  }
  return {mu, theta, counts, fitted, sweeps, converged};
}

/** Pearson correlation weighted by each level's support, so a value carrying five hundred compounds
 *  is not outvoted by one carrying three. */
function weightedCorrelation(xs: number[], ys: number[], weights: number[]): number | null {
  let total = 0;
  let mx = 0;
  let my = 0;
  for (let i = 0; i < xs.length; i++) {
    total += weights[i];
    mx += weights[i] * xs[i];
    my += weights[i] * ys[i];
  }
  mx /= total;
  my /= total;
  let cov = 0;
  let vx = 0;
  let vy = 0;
  for (let i = 0; i < xs.length; i++) {
    const dx = xs[i] - mx;
    const dy = ys[i] - my;
    cov += weights[i] * dx * dy;
    vx += weights[i] * dx * dx;
    vy += weights[i] * dy * dy;
  }
  return (vx === 0 || vy === 0) ? null : cov / Math.sqrt(vx * vy);
}

export function fitRoleEffects(design: RoleDesign): RoleFit | null {
  const roleCount = design.names.length;
  if (roleCount < 2)
    return null;
  const n = design.activity.length;

  const codes: Int32Array[] = [];
  const levelValues: string[][] = [];
  for (let r = 0; r < roleCount; r++) {
    const index = new Map<string, number>();
    const coded = new Int32Array(n);
    const values: string[] = [];
    for (let k = 0; k < n; k++) {
      let code = index.get(design.values[r][k]);
      if (code === undefined) {
        code = values.length;
        index.set(design.values[r][k], code);
        values.push(design.values[r][k]);
      }
      coded[k] = code;
    }
    codes.push(coded);
    levelValues.push(values);
  }

  // Offsets are identified only within a connected component of the design: two blocks sharing no
  // value in any role fix their offsets against unrelated baselines, so putting one block's offset
  // beside the other's is arithmetic over incomparable quantities.
  const base: number[] = [];
  let nodes = 0;
  for (let r = 0; r < roleCount; r++) {
    base.push(nodes);
    nodes += levelValues[r].length;
  }
  const parent = new Int32Array(nodes).map((_v, i) => i);
  const find = (x: number): number => {
    while (parent[x] !== x)
      x = parent[x] = parent[parent[x]];
    return x;
  };
  for (let k = 0; k < n; k++) {
    for (let r = 1; r < roleCount; r++) {
      const a = find(codes[0][k]);
      const b = find(base[r] + codes[r][k]);
      if (a !== b)
        parent[a] = b;
    }
  }
  const perComponent = new Map<number, number>();
  for (let k = 0; k < n; k++) {
    const root = find(codes[0][k]);
    perComponent.set(root, (perComponent.get(root) ?? 0) + 1);
  }
  let largest = -1;
  let largestSize = 0;
  for (const [root, size] of perComponent) {
    if (size > largestSize) {
      largestSize = size;
      largest = root;
    }
  }
  const keep: number[] = [];
  for (let k = 0; k < n; k++) {
    if (find(codes[0][k]) === largest)
      keep.push(k);
  }

  // Every fit runs over the kept observations only, so a level living only in a dropped block keeps a
  // count of zero and stays out of every sum below rather than sitting at a coefficient nothing measured.
  const m = keep.length;
  const levelCounts = levelValues.map((values) => values.length);
  const y = Float64Array.from(design.activity);
  const mol = design.molIdx;
  const full = backfit(codes, levelCounts, y, Int32Array.from(keep));
  const ranked = full.counts.map((counts) =>
    counts.reduce((total, support) => total + (support >= design.minSupport ? 1 : 0), 0));
  if (!ranked.some((count) => count >= 2))
    return null;

  let sse = 0;
  for (const k of keep)
    sse += (y[k] - full.fitted[k]) ** 2;
  // Effective rather than nominal degrees of freedom: the prior costs a level less than a whole
  // parameter, and the difference grows with the number of thinly supported levels.
  let edf = 1;
  for (let r = 0; r < roleCount; r++) {
    for (const support of full.counts[r])
      edf += support / (support + ROLE_FIT_PRIOR);
  }
  const residualSd = Math.sqrt(sse / Math.max(1, m - edf));

  // Folds key on the source row, never on position: fragments arrive in worker-completion order, so
  // a fold keyed on the walk index would give a different R² for identical input.
  let ssRes = 0;
  for (let fold = 0; fold < ROLE_CV_FOLDS; fold++) {
    const train: number[] = [];
    const test: number[] = [];
    for (const k of keep)
      (mol[k] % ROLE_CV_FOLDS === fold ? test : train).push(k);
    if (test.length === 0)
      continue;
    const trained = backfit(codes, levelCounts, y, Int32Array.from(train));
    for (const i of test) {
      let predicted = trained.mu;
      // A level the training folds never saw is left at zero, which by the centring is that role's
      // library average — so the rest of the compound still informs the prediction.
      for (let r = 0; r < roleCount; r++)
        predicted += trained.theta[r][codes[r][i]];
      ssRes += (y[i] - predicted) ** 2;
    }
  }
  let cvR2: number | null = null;
  let cvRmse: number | null = null;
  if (m >= ROLE_MIN_CV_POINTS) {
    const cvMean = keep.reduce((sum, k) => sum + y[k], 0) / m;
    const ssTot = keep.reduce((sum, k) => sum + (y[k] - cvMean) ** 2, 0);
    if (ssTot !== 0) {
      cvR2 = 1 - ssRes / ssTot;
      cvRmse = Math.sqrt(ssRes / m);
    }
  }

  // Cross-validated R² checks prediction, and prediction is nearly unique even where the split of
  // credit between roles is not. Refitting each half of the table is the only check on that split.
  const halves = [0, 1].map((half) =>
    backfit(codes, levelCounts, y, Int32Array.from(keep.filter((k) => mol[k] % 2 === half))));

  const error = cvRmse ?? residualSd;
  const dir = design.higherIsBetter ? 1 : -1;
  const roles: RoleSummary[] = [];
  for (let r = 0; r < roleCount; r++) {
    const counts = full.counts[r];
    const theta = full.theta[r];
    const observed = counts.filter((count) => count > 0).length;
    const levels: RoleLevel[] = [];
    let rawVar = 0;
    const xs: number[] = [];
    const ys: number[] = [];
    const weights: number[] = [];
    for (let v = 0; v < counts.length; v++) {
      rawVar += counts[v] * theta[v] * theta[v];
      if (counts[v] >= design.minSupport)
        levels.push({value: levelValues[r][v], coef: theta[v], n: counts[v]});
      if (halves[0].counts[r][v] >= design.minSupport && halves[1].counts[r][v] >= design.minSupport) {
        xs.push(halves[0].theta[r][v]);
        ys.push(halves[1].theta[r][v]);
        weights.push(counts[v]);
      }
    }
    levels.sort((a, b) => dir * (b.coef - a.coef));
    rawVar /= m;
    // The raw count-weighted sd carries an estimation-noise term σ²·levels/observations that is
    // negligible for a twelve-level role and dominant for a six-hundred-level one, so subtracting it
    // is what makes two roles of one table comparable at all.
    const spread = Math.sqrt(Math.max(0, rawVar - residualSd ** 2 * observed / m));

    let cut = 0;
    let gap = 0;
    for (let i = 1; i < levels.length; i++) {
      const step = dir * (levels[i - 1].coef - levels[i].coef);
      if (step > gap) {
        gap = step;
        cut = i;
      }
    }
    let split: RoleSplit | null = null;
    // A lone outlier at either end is one value, not a group, so the widest gap is reported only
    // where it partitions the role rather than clipping it.
    if (gap > error && cut >= 2 && levels.length - cut >= 2) {
      const side = (from: number, to: number): {mean: number, total: number} => {
        let weighted = 0;
        let total = 0;
        for (let i = from; i < to; i++) {
          weighted += levels[i].n * levels[i].coef;
          total += levels[i].n;
        }
        return {mean: weighted / total, total};
      };
      const hi = side(0, cut);
      const lo = side(cut, levels.length);
      split = {gap, hiCount: cut, hiMean: hi.mean, hiN: hi.total,
        loCount: levels.length - cut, loMean: lo.mean, loN: lo.total};
    }

    // Two comparable levels correlate at exactly ±1 whatever the halves hold, so a two-value role
    // would either wear the strongest possible repeatability chip or be refused outright, on nothing.
    roles.push({name: design.names[r], levels, fitted: observed,
      thin: observed - levels.length, spread,
      repeat: xs.length < 3 ? null : weightedCorrelation(xs, ys, weights), split});
  }
  roles.sort((a, b) => b.spread - a.spread);

  return {roles, mean: full.mu, compounds: m, dropped: n - m, sweeps: full.sweeps,
    converged: full.converged, residualSd, cvR2, cvRmse};
}
