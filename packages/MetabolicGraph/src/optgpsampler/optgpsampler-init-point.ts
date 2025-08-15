/*
  OptGPSampler (TypeScript)
  ---------------------------------
  A simplified, dependency‑free implementation inspired by COBRApy's OptGPSampler,
  now with an **automatic feasible point finder** (so you can omit x0).

  Notes
  -----
  • Targets the same *idea* as OptGP/ACHR: centering hit‑and‑run with short chains.
  • Works for polytopes defined by:
      A x <= b
      Aeq x = beq
      lb <= x <= ub
    and **can compute a feasible starting point** if not provided.
  • Equality feasibility is preserved by projecting directions to the nullspace of Aeq.
  • Sampling along a line uses a uniform step within feasible [alpha_min, alpha_max].
  • "nproj" controls how many internal steps are taken for each returned sample.
  • "thinning" controls how many steps are skipped between recorded chain states.
  • RNG is seedable and deterministic.

  Feasible point finder
  ---------------------
  We implement a lightweight, dependency‑free Phase‑I style search:
    1) Compute a particular solution x_p for Aeq x = beq (minimum‑norm projection); if Aeq is absent, start from mid‑bounds.
    2) Alternating projection between equality manifold and box bounds.
    3) Reduce inequality violations A x <= b by projected gradient steps (in the Aeq nullspace), with backtracking and box clipping.
  If the region is feasible and reasonably conditioned, this converges quickly in practice.
  (It is not a full LP solver, but behaves similarly to COBRApy's practice of finding a feasible LP solution.)
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
  next(): number { return this.nextUint32() / 0x100000000; } // [0,1)
  normal(): number {
    // Box-Muller
    let u = 0, v = 0;
    while (u === 0) u = this.next();
    while (v === 0) v = this.next();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  }
}

// ----------------------------- Utilities: LA ------------------------------
function dot(a: number[], b: number[]): number { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }
function add(a: number[], b: number[]): number[] { const n = a.length; const out = new Array(n); for (let i = 0; i < n; i++) out[i] = a[i] + b[i]; return out; }
function sub(a: number[], b: number[]): number[] { const n = a.length; const out = new Array(n); for (let i = 0; i < n; i++) out[i] = a[i] - b[i]; return out; }
function scale(a: number[], c: number): number[] { const n = a.length; const out = new Array(n); for (let i = 0; i < n; i++) out[i] = a[i] * c; return out; }
function clone(a: number[]): number[] { return a.slice(); }
function norm2(a: number[]): number { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * a[i]; return Math.sqrt(s); }
function normInf(a: number[]): number { let m = 0; for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i])); return m; }
function clamp(x: number[], lo: number[], hi: number[]): number[] { const n = x.length; const y = new Array(n); for (let i=0;i<n;i++){ y[i] = Math.min(hi[i], Math.max(lo[i], x[i])); } return y; }

function matVec(A: number[][], x: number[]): number[] {
  const m = A.length, n = x.length; const out = new Array(m);
  for (let i = 0; i < m; i++) { let s = 0; const row = A[i]; for (let j = 0; j < n; j++) s += row[j] * x[j]; out[i] = s; }
  return out;
}
function transpose(A: number[][]): number[][] {
  const m = A.length, n = A[0]?.length ?? 0;
  const AT: number[][] = Array.from({ length: n }, () => new Array(m).fill(0));
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) AT[j][i] = A[i][j];
  return AT;
}

// Solve linear system via Gauss-Jordan (simple, no deps)
function solveLinear(Ain: number[][], b: number[]): number[] {
  const n = b.length;
  const A = Ain.map(r => r.slice());
  const x = b.slice();
  if (A.length !== n) throw new Error('solveLinear expects square matrix');
  // Augment
  for (let i = 0; i < n; i++) A[i].push(x[i]);
  // Elimination
  for (let i = 0; i < n; i++) {
    // pivot
    let piv = i, best = Math.abs(A[i][i]);
    for (let r = i + 1; r < n; r++) { const v = Math.abs(A[r][i]); if (v > best) { best = v; piv = r; } }
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
  return sub(v, ATL);
}

// Compute a particular solution of Aeq x = beq with minimum 2-norm (via normal equations)
function eqParticular(Aeq: number[][], beq: number[]): number[] {
  const AT = transpose(Aeq); // n x m
  // Solve (Aeq^T Aeq) x = Aeq^T beq
  const ATA: number[][] = Array.from({ length: AT.length }, () => new Array(AT.length).fill(0));
  for (let i = 0; i < AT.length; i++) {
    for (let j = 0; j < AT.length; j++) {
      let s = 0; for (let k = 0; k < Aeq.length; k++) s += AT[i][k] * Aeq[k][j];
      ATA[i][j] = s;
    }
  }
  const ATb = new Array(AT.length).fill(0);
  for (let i = 0; i < AT.length; i++) { let s = 0; for (let k = 0; k < Aeq.length; k++) s += AT[i][k] * beq[k]; ATb[i] = s; }
  return solveLinear(ATA, ATb);
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
  // A feasible starting point (OPTIONAL now)
  x0?: number[];
}

export interface OptGPOptions {
  thinning?: number;   // return every k steps (default 100)
  nproj?: number;      // short chain length per sample (default 8)
  warmupSteps?: number; // steps to estimate center (default 5_000)
  seed?: number;       // RNG seed
  // Feasible point finder controls
  feasTol?: number;            // 1e-7
  feasMaxIter?: number;        // 5000
  feasIneqWeight?: number;     // gradient weight for inequalities (1.0)
  feasStep?: number;           // initial step size (1.0)
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
  private feasTol: number;
  private feasMaxIter: number;
  private feasIneqWeight: number;
  private feasStep: number;

  constructor(model: LinearPolytope, opts: OptGPOptions = {}) {
    // infer n from bounds or Aeq or A or x0
    const nCand = model.x0?.length ?? model.lb?.length ?? model.ub?.length ?? model.Aeq?.[0]?.length ?? model.A?.[0]?.length;
    if (!nCand || nCand <= 0) throw new Error('Cannot infer variable dimension n; provide lb/ub/Aeq/A or x0');
    this.n = nCand;

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
    this.feasTol = opts.feasTol ?? 1e-7;
    this.feasMaxIter = Math.max(1, Math.floor(opts.feasMaxIter ?? 5000));
    this.feasIneqWeight = opts.feasIneqWeight ?? 1.0;
    this.feasStep = opts.feasStep ?? 1.0;

    // Dimension checks
    const checkVec = (v: number[] | undefined | null, name: string) => { if (!v) return; if (v.length !== this.n) throw new Error(`${name} has wrong length`); };
    checkVec(model.lb, 'lb'); checkVec(model.ub, 'ub'); checkVec(model.x0, 'x0');
    if (this.A && this.b) {
      if (this.A.length !== this.b.length) throw new Error('A and b dimension mismatch');
      this.A.forEach((row, i) => { if (row.length !== this.n) throw new Error(`A[${i}] length mismatch`); });
    }
    if (this.Aeq && this.beq) {
      if (this.Aeq.length !== this.beq.length) throw new Error('Aeq and beq dimension mismatch');
      this.Aeq.forEach((row, i) => { if (row.length !== this.n) throw new Error(`Aeq[${i}] length mismatch`); });
    }

    // Initialize x: use given x0 if feasible, else find one
    const initial = model.x0?.slice();
    if (initial && this.isFeasible(initial)) {
      this.x = initial;
    } else {
      this.x = this.findFeasiblePoint(initial ?? null);
    }
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

  // ---------------- Feasible point finder ----------------
  private findFeasiblePoint(seed: number[] | null): number[] {
    // Start from mid-bounds if available, else zeros
    let x = new Array(this.n).fill(0);
    const mid = new Array(this.n);
    for (let i = 0; i < this.n; i++) {
      const lo = this.lb[i], hi = this.ub[i];
      mid[i] = (isFinite(lo) && isFinite(hi)) ? 0.5 * (lo + hi) : 0;
    }
    x = seed ?? mid;

    // If equality constraints: compute particular solution and start from there
    if (this.Aeq && this.Aeq.length > 0) {
      try {
        const xp = eqParticular(this.Aeq, this.beq as number[]);
        x = xp;
      } catch { /* fall back to mid */ }
    }

    // Alternating projection between equality manifold and bounds
    const projectEqualities = (vec: number[]) => {
      if (!this.Aeq || this.Aeq.length === 0) return vec;
      const r = sub(matVec(this.Aeq, vec), this.beq as number[]); // residual
      const corr = projectToNullspace(this.Aeq, r); // A^T(A A^T)^{-1} (A x - b) in column space, but we want to subtract it
      // Above helper returns v - A^T (A A^T)^{-1} (A v). We need to *fix* vec: vec' = vec - A^T(...) r
      // Construct A^T(...)r directly:
      const m = this.Aeq.length; const AT = transpose(this.Aeq);
      // Solve (A A^T) y = r, then AT y
      const AAT: number[][] = Array.from({ length: m }, (_, i) => {
        const row: number[] = new Array(m);
        for (let j = 0; j < m; j++) { let s = 0; for (let k = 0; k < this.n; k++) s += this.Aeq![i][k] * this.Aeq![j][k]; row[j] = s; }
        return row;
      });
      let y: number[];
      try { y = solveLinear(AAT, r); } catch { return vec; }
      const ATy = new Array(this.n).fill(0);
      for (let i = 0; i < this.n; i++) { let s = 0; for (let j = 0; j < m; j++) s += AT[i][j] * y[j]; ATy[i] = s; }
      return sub(vec, ATy);
    };

    // Iterate alternating projections
    for (let it = 0; it < Math.min(200, this.feasMaxIter); it++) {
      x = clamp(x, this.lb, this.ub);
      const xOld = x;
      x = projectEqualities(x);
      if (normInf(sub(x, xOld)) < this.feasTol) break;
    }

    // Reduce inequality violations with projected gradient descent in nullspace of Aeq
    let step = this.feasStep;
    for (let it = 0; it < this.feasMaxIter; it++) {
      // violation vector v >= 0 where A x - b > 0
      let violated = false;
      let grad = new Array(this.n).fill(0);
      if (this.A) {
        const Ax = matVec(this.A, x);
        for (let i = 0; i < Ax.length; i++) {
          const val = Ax[i] - (this.b as number[])[i];
          if (val > 0) {
            violated = true;
            const row = this.A[i];
            for (let j = 0; j < this.n; j++) grad[j] += val * row[j];
          }
        }
      }
      if (!violated) {
        // final tighten equality + bounds
        x = clamp(projectEqualities(x), this.lb, this.ub);
        if (this.isFeasible(x)) return x;
      }
      // Move in negative gradient, keep within equality nullspace
      if (norm2(grad) < 1e-16) break;
      let d = scale(grad, -this.feasIneqWeight);
      d = projectToNullspace(this.Aeq ?? null, d);

      // backtracking line search with clipping
      let accepted = false;
      for (let bt = 0; bt < 20; bt++) {
        const cand = clamp(add(x, scale(d, step)), this.lb, this.ub);
        const better = this.violation(cand) < this.violation(x) - 1e-12;
        if (better) { x = projectEqualities(cand); accepted = true; break; }
        step *= 0.5;
      }
      if (!accepted) break;
      // check small overall violation
      if (this.violation(x) < this.feasTol) {
        x = clamp(projectEqualities(x), this.lb, this.ub);
        if (this.isFeasible(x)) return x;
      }
    }

    if (!this.isFeasible(x)) {
      throw new Error('Failed to find a feasible starting point. Consider providing x0 or loosening constraints.');
    }
    return x;
  }

  private violation(x: number[]): number {
    let v = 0;
    // bounds
    for (let i = 0; i < this.n; i++) {
      if (x[i] < this.lb[i]) v += (this.lb[i] - x[i]);
      if (x[i] > this.ub[i]) v += (x[i] - this.ub[i]);
    }
    // inequalities
    if (this.A) {
      const Ax = matVec(this.A, x);
      for (let i = 0; i < Ax.length; i++) v += Math.max(0, Ax[i] - (this.b as number[])[i]);
    }
    // equalities
    if (this.Aeq) {
      const Ax = matVec(this.Aeq, x);
      for (let i = 0; i < Ax.length; i++) v += Math.abs(Ax[i] - (this.beq as number[])[i]);
    }
    return v;
  }

  // ---------------- Warmup + Sampling ----------------
  private warmup(): void {
    if (this.center) return;
    let x = this.x.slice();
    let mean = new Array(this.n).fill(0);
    let cnt = 0;
    for (let t = 0; t < this.warmupSteps; t++) {
      x = this.singleStep(x);
      if ((t + 1) % this.thinning === 0) { cnt++; for (let i = 0; i < this.n; i++) mean[i] += (x[i] - mean[i]) / cnt; }
    }
    this.center = cnt > 0 ? mean : x.slice();
    this.x = x;
  }

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
        if (Math.abs(c) < 1e-14) { if (rhs < -1e-12) return { amin: 1, amax: 0 }; continue; }
        const a = rhs / c;
        if (c > 0) { if (a < amax) amax = a; } else { if (a > amin) amin = a; }
      }
    }
    return { amin, amax };
  }

  private singleStep(x: number[]): number[] {
    const c = this.center ?? x; // before warmup, center is unknown
    const g = new Array(this.n);
    for (let i = 0; i < this.n; i++) g[i] = this.rng.normal();
    let d = add(sub(c, x), g);
    d = projectToNullspace(this.Aeq ?? null, d);

    const { amin, amax } = this.stepInterval(x, d);
    if (!(amin < amax) || !isFinite(amin) || !isFinite(amax)) {
      // Fallback to a pure random nullspace direction
      let r = new Array(this.n); for (let i = 0; i < this.n; i++) r[i] = this.rng.normal();
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

    const takeStep = () => { x = this.singleStep(x); };

    for (let s = 0; s < N; s++) {
      for (let j = 0; j < this.nproj; j++) {
        for (let k = 0; k < this.thinning; k++) takeStep();
      }
      out.push(clone(x));
    }
    this.x = x;
    return out;
  }
}

// ----------------------------- Example usage ------------------------------

// Suppose you have Sv = 0, with bounds lb <= v <= ub (COBRA-style)
/*const S = [ [-2, 1, 0], [0, 1, -2] ];
const beq = [0, 0];
const lb = [0, 0, 0];
const ub = [10, 10, 10];

// You can now omit x0 — the sampler will find one automatically
const sampler = new OptGPSampler({ Aeq: S, beq, lb, ub }, {
  thinning: 50,
  nproj: 8,
  warmupSteps: 2000,
  seed: 42,
  feasTol: 1e-8,
  feasMaxIter: 4000,
});

const samples = sampler.sample(1000);
console.log(samples.length, 'samples, first:', samples[220]);*/

const S = [ [-2, 1, 0, 1], [0, 1, -2, 1] ];
const beq = [0, 0];
const lb = [-10, -10, -10, -10];
const ub = [10, 10, 10, 10];

// You can now omit x0 — the sampler will find one automatically
const sampler = new OptGPSampler({ Aeq: S, beq, lb, ub }, {
  thinning: 50,
  nproj: 8,
  warmupSteps: 2000,
  seed: 42,
  feasTol: 1e-8,
  feasMaxIter: 4000,
});

const samples = sampler.sample(1000);
console.log(samples.length, 'samples, first:', samples[220]);

