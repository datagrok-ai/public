/*
  OptGPSampler (TypeScript)
  ---------------------------------
  A simplified, dependency‑free implementation inspired by COBRApy's OptGPSampler.

  Notes
  -----
  • Targets the same *idea* as OptGP/ACHR: center‑ing hit‑and‑run with short chains.
  • Works for polytopes defined by:
      A x <= b
      Aeq x = beq
      lb <= x <= ub
    and requires a feasible starting point x0.
  • Equality feasibility is preserved by projecting directions to the nullspace of Aeq.
  • Sampling along a line uses a uniform step within feasible [alpha_min, alpha_max].
  • "nproj" controls how many internal steps are taken for each returned sample.
  • "thinning" controls how many steps are skipped between recorded chain states.
  • RNG is seedable and deterministic.

  This is not a drop‑in reimplementation of COBRApy; it’s a practical TypeScript sampler
  with similar behavior and API surface for Node or browser builds.
*/

// ----------------------------- Utilities: RNG -----------------------------
class XorShift32 {
  private state: number;
  constructor(seed = 123456789) {
    if (seed === 0) seed = 1;
    this.state = seed >>> 0;
  }
  nextUint32(): number {
    // xorshift32
    let x = this.state;
    x ^= x << 13; x >>>= 0;
    x ^= x >>> 17; x >>>= 0;
    x ^= x << 5;  x >>>= 0;
    this.state = x >>> 0;
    return this.state;
  }
  next(): number {
    // [0,1)
    return this.nextUint32() / 0x100000000;
  }
  normal(): number {
    // Box-Muller
    let u = 0, v = 0;
    while (u === 0) u = this.next();
    while (v === 0) v = this.next();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  }
}

// ----------------------------- Utilities: LA ------------------------------
function dot(a: number[], b: number[]): number {
  let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s;
}
function add(a: number[], b: number[]): number[] {
  const n = a.length; const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a[i] + b[i];
  return out;
}
function sub(a: number[], b: number[]): number[] {
  const n = a.length; const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a[i] - b[i];
  return out;
}
function scale(a: number[], c: number): number[] {
  const n = a.length; const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a[i] * c;
  return out;
}
function clone(a: number[]): number[] { return a.slice(); }

function matVec(A: number[][], x: number[]): number[] {
  const m = A.length, n = x.length; const out = new Array(m);
  for (let i = 0; i < m; i++) {
    let s = 0; const row = A[i];
    for (let j = 0; j < n; j++) s += row[j] * x[j];
    out[i] = s;
  }
  return out;
}
function transpose(A: number[][]): number[][] {
  const m = A.length, n = A[0]?.length ?? 0;
  const AT: number[][] = Array.from({ length: n }, () => new Array(m).fill(0));
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) AT[j][i] = A[i][j];
  return AT;
}

// Solve symmetric positive definite or well-conditioned system via Gauss-Jordan
function solveLinear(Ain: number[][], b: number[]): number[] {
  const n = b.length;
  const A = Ain.map(r => r.slice());
  const x = b.slice();
  // Augment
  for (let i = 0; i < n; i++) A[i].push(x[i]);
  // Elimination
  for (let i = 0; i < n; i++) {
    // pivot
    let piv = i, best = Math.abs(A[i][i]);
    for (let r = i + 1; r < n; r++) {
      const v = Math.abs(A[r][i]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < 1e-12) throw new Error('Singular matrix in solveLinear');
    if (piv !== i) { const tmp = A[i]; A[i] = A[piv]; A[piv] = tmp; }
    // normalize
    const diag = A[i][i];
    for (let c = i; c <= n; c++) A[i][c] /= diag;
    // eliminate
    for (let r = 0; r < n; r++) if (r !== i) {
      const f = A[r][i];
      for (let c = i; c <= n; c++) A[r][c] -= f * A[i][c];
    }
  }
  // Extract solution
  const sol = new Array(n);
  for (let i = 0; i < n; i++) sol[i] = A[i][n];
  return sol;
}

// Project vector v to the nullspace of Aeq: v' = v - Aeq^T (Aeq Aeq^T)^-1 (Aeq v)
function projectToNullspace(Aeq: number[][] | null, v: number[]): number[] {
  if (!Aeq || Aeq.length === 0) return v.slice();
  const m = Aeq.length; const Av = matVec(Aeq, v); // m
  const AAT: number[][] = Array.from({ length: m }, (_, i) => {
    const row: number[] = new Array(m);
    for (let j = 0; j < m; j++) {
      let s = 0; const ai = Aeq[i], aj = Aeq[j];
      for (let k = 0; k < ai.length; k++) s += ai[k] * aj[k];
      row[j] = s;
    }
    return row;
  });
  const lambda = solveLinear(AAT, Av); // m
  const AT = transpose(Aeq);
  // AT * lambda
  const ATL = new Array(AT.length).fill(0);
  for (let i = 0; i < AT.length; i++) {
    let s = 0; const col = AT[i];
    for (let j = 0; j < m; j++) s += col[j] * lambda[j];
    ATL[i] = s;
  }
  const proj = sub(v, ATL);
  return proj;
}

// ----------------------------- Model types ------------------------------
export interface LinearPolytope {
  // Inequalities: A x <= b
  A?: number[][]; // (mi x n)
  b?: number[];   // (mi)
  // Equalities: Aeq x = beq
  Aeq?: number[][]; // (me x n)
  beq?: number[];   // (me)
  // Bounds: lb <= x <= ub (componentwise)
  lb?: number[];
  ub?: number[];
  // A feasible starting point (REQUIRED)
  x0: number[];
}

export interface OptGPOptions {
  thinning?: number;   // return every k steps (default 100)
  nproj?: number;      // short chain length per sample (default 8)
  warmupSteps?: number; // steps to estimate center (default 5_000)
  seed?: number;       // RNG seed
}

// ----------------------------- Sampler ------------------------------
export class OptGPSampler {
  private readonly n: number;
  private readonly A: number[][] | null;
  private readonly b: number[] | null;
  private readonly Aeq: number[][] | null;
  private readonly beq: number[] | null;
  private readonly lb: number[];
  private readonly ub: number[];
  private readonly rng: XorShift32;
  private center: number[] | null = null;
  private thinning: number;
  private nproj: number;
  private warmupSteps: number;
  private x: number[]; // current state

  constructor(model: LinearPolytope, opts: OptGPOptions = {}) {
    this.n = model.x0.length;
    this.A = model.A ?? null;
    this.b = model.b ?? null;
    this.Aeq = model.Aeq ?? null;
    this.beq = model.beq ?? null;
    this.lb = model.lb ?? new Array(this.n).fill(-Infinity);
    this.ub = model.ub ?? new Array(this.n).fill(+Infinity);
    this.rng = new XorShift32(opts.seed ?? 1234);
    this.thinning = Math.max(1, Math.floor(opts.thinning ?? 100));
    this.nproj = Math.max(1, Math.floor(opts.nproj ?? 8));
    this.warmupSteps = Math.max(0, Math.floor(opts.warmupSteps ?? 5000));

    // Verify dimensions
    const checkVec = (v: number[] | undefined | null, name: string) => {
      if (!v) return;
      if (v.length !== this.n) throw new Error(`${name} has wrong length`);
    };
    checkVec(model.lb, 'lb'); checkVec(model.ub, 'ub');
    if (this.A && this.b) {
      if (this.A.length !== this.b.length) throw new Error('A and b dimension mismatch');
      this.A.forEach((row, i) => { if (row.length !== this.n) throw new Error(`A[${i}] length mismatch`); });
    }
    if (this.Aeq && this.beq) {
      if (this.Aeq.length !== this.beq.length) throw new Error('Aeq and beq dimension mismatch');
      this.Aeq.forEach((row, i) => { if (row.length !== this.n) throw new Error(`Aeq[${i}] length mismatch`); });
    }

    // Check that x0 is feasible
    this.x = model.x0.slice();
    if (!this.isFeasible(this.x)) throw new Error('x0 is not feasible for the given constraints.');
  }

  // Check feasibility
  private isFeasible(x: number[]): boolean {
    for (let i = 0; i < this.n; i++) {
      if (x[i] < this.lb[i] - 1e-8 || x[i] > this.ub[i] + 1e-8) return false;
    }
    if (this.A) {
      const Ax = matVec(this.A, x);
      for (let i = 0; i < Ax.length; i++) if (Ax[i] - (this.b as number[])[i] > 1e-8) return false;
    }
    if (this.Aeq) {
      const Ax = matVec(this.Aeq, x);
      for (let i = 0; i < Ax.length; i++) if (Math.abs(Ax[i] - (this.beq as number[])[i]) > 1e-8) return false;
    }
    return true;
  }

  // Warmup to estimate the center using hit-and-run with nullspace projection
  private warmup(): void {
    if (this.center) return;
    let x = this.x.slice();
    let mean = new Array(this.n).fill(0);
    let cnt = 0;
    for (let t = 0; t < this.warmupSteps; t++) {
      x = this.singleStep(x);
      // online mean after thinning to avoid heavy autocorrelation
      if ((t + 1) % this.thinning === 0) {
        cnt++;
        for (let i = 0; i < this.n; i++) mean[i] += (x[i] - mean[i]) / cnt;
      }
    }
    // Fallback if warmupSteps == 0
    this.center = cnt > 0 ? mean : x.slice();
    this.x = x;
  }

  // Compute feasible step interval along direction d from point x
  private stepInterval(x: number[], d: number[]): { amin: number; amax: number } {
    let amin = -Infinity, amax = +Infinity;
    // Bounds
    for (let i = 0; i < this.n; i++) {
      if (Math.abs(d[i]) < 1e-14) continue;
      const a1 = (this.lb[i] - x[i]) / d[i];
      const a2 = (this.ub[i] - x[i]) / d[i];
      const lo = Math.min(a1, a2), hi = Math.max(a1, a2);
      if (lo > amin) amin = lo;
      if (hi < amax) amax = hi;
    }
    // Inequalities: A x <= b  =>  a * (A d) <= b - A x
    if (this.A) {
      const Ad = matVec(this.A, d);
      const Ax = matVec(this.A, x);
      for (let i = 0; i < Ad.length; i++) {
        const c = Ad[i];
        const rhs = (this.b as number[])[i] - Ax[i];
        if (Math.abs(c) < 1e-14) {
          if (rhs < -1e-12) return { amin: 1, amax: 0 }; // infeasible line
          continue;
        }
        const a = rhs / c;
        if (c > 0) { // a <= rhs/c
          if (a < amax) amax = a;
        } else { // a >= rhs/c
          if (a > amin) amin = a;
        }
      }
    }
    return { amin, amax };
  }

  // One projected hit-and-run step with centering perturbation
  private singleStep(x: number[]): number[] {
    const c = this.center ?? x; // before warmup, center is unknown
    // direction: (c - x) + random normal vector
    const g = new Array(this.n);
    for (let i = 0; i < this.n; i++) g[i] = this.rng.normal();
    let d = add(sub(c, x), g);
    d = projectToNullspace(this.Aeq ?? null, d);

    // If equality constraints exist, also ensure x stays on Aeq x = beq by stepping in nullspace only
    // Compute feasible interval and sample uniformly
    const { amin, amax } = this.stepInterval(x, d);
    if (!(amin < amax) || !isFinite(amin) || !isFinite(amax)) {
      // Degenerate direction; try a pure random nullspace direction
      let r = new Array(this.n);
      for (let i = 0; i < this.n; i++) r[i] = this.rng.normal();
      d = projectToNullspace(this.Aeq ?? null, r);
      const seg = this.stepInterval(x, d);
      if (!(seg.amin < seg.amax)) return x; // give up, rare
      const a = seg.amin + this.rng.next() * (seg.amax - seg.amin);
      return add(x, scale(d, a));
    }
    const a = amin + this.rng.next() * (amax - amin);
    return add(x, scale(d, a));
  }

  // Public: generate N samples
  sample(N: number): number[][] {
    if (N <= 0) return [];
    this.warmup();
    const out: number[][] = [];
    let x = this.x.slice();
    let stepsSinceRecord = 0;

    const takeStep = () => { x = this.singleStep(x); stepsSinceRecord++; };

    for (let s = 0; s < N; s++) {
      // short internal chain (nproj steps), record only last point
      for (let j = 0; j < this.nproj; j++) {
        // maintain thinning across steps so correlation is reduced
        for (let k = 0; k < this.thinning; k++) takeStep();
      }
      out.push(clone(x));
      stepsSinceRecord = 0;
    }
    // update internal state
    this.x = x;
    return out;
  }
}

// ----------------------------- Example usage ------------------------------

// Suppose you have Sv = 0, with bounds lb <= v <= ub
const S = [ [2, -1, 0], [0, 1, -2] ];
const beq = [0, 0];
const lb = [0, 0, 0];
const ub = [10, 10, 10];
const x0 = [1, 2, 1]; // feasible (S x0 = 0)

const sampler = new OptGPSampler({ Aeq: S, beq, lb, ub, x0 }, {
  thinning: 50,
  nproj: 8,
  warmupSteps: 2000,
  seed: 42,
});

const samples = sampler.sample(1000);
/*console.log(samples.length, 'samples, first:', samples[0]);
console.log(samples.length, 'samples, first:', samples[1]);
console.log(samples.length, 'samples, first:', samples[2]);*/

samples.forEach((s) => console.log(s));
