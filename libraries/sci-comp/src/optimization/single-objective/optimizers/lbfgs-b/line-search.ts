/**
 * Moré–Thuente line search enforcing the strong Wolfe conditions.
 *
 * Port of MINPACK-2 `dcsrch` + `dcstep` (Moré & Thuente 1994). The
 * line search seeks a step α > 0 along the one-dimensional function
 *
 *    φ(α) = f(x + α·d),    φ'(α) = ∇f(x + α·d)ᵀ d,
 *
 * satisfying both the sufficient-decrease (Armijo) and curvature
 * conditions:
 *
 *    φ(α) ≤ φ(0) + c₁·α·φ'(0)     (ftol)
 *    |φ'(α)| ≤ c₂·|φ'(0)|          (gtol)
 *
 * Internally, stage 1 works with the auxiliary
 *    ψ(α) = φ(α) − φ(0) − c₁·α·φ'(0),
 * switching to φ once a low-enough function value is seen and the slope
 * is non-negative. A safeguarded cubic/quadratic step (`dcstep`) handles
 * the four geometric cases that arise while bracketing the minimum.
 *
 * Reference:
 *  - J. J. Moré, D. J. Thuente (1994), *Line Search Algorithms with
 *    Guaranteed Sufficient Decrease*, ACM TOMS 20(3):286–307.
 *  - MINPACK-2 source: `dcsrch.f`, `dcstep.f`.
 */

/* ================================================================== */
/*  Public types                                                       */
/* ================================================================== */

/** Evaluator for φ(α) = f(x + αd) and φ'(α) = ∇f(x + αd)ᵀd. */
export type PhiEval = (alpha: number) => {phi: number; phiPrime: number};

/** Async variant of {@link PhiEval}. */
export type PhiEvalAsync = (alpha: number) => Promise<{phi: number; phiPrime: number}>;

/** Line-search tuning knobs. */
export interface LineSearchParams {
  /** Armijo sufficient-decrease constant c₁ ∈ (0, 1). */
  ftol: number;
  /** Curvature constant c₂ ∈ (c₁, 1). */
  gtol: number;
  /** Relative interval-width tolerance. */
  xtol: number;
  /** Maximum number of φ/φ' evaluations. */
  maxSteps: number;
  /** Upper cap on α (typically `maxFeasibleStep`). */
  stpMax: number;
  /** Lower cap on α. Default 1e-20 (MINPACK-2). */
  stpMin?: number;
}

/** Outcome classification of the line search. */
export type LineSearchStatus =
  | 'converged' // strong Wolfe satisfied
  | 'warning_rounding' // rounding errors prevent further progress
  | 'warning_xtol' // interval width below xtol·stmax
  | 'warning_stpmax' // stp hit stpMax with descent still valid
  | 'warning_stpmin' // stp hit stpMin with descent still valid
  | 'max_steps_reached' // exhausted maxSteps
  | 'error'; // invalid inputs (defensive)

export interface LineSearchResult {
  status: LineSearchStatus;
  /** Accepted (or best-so-far) step. */
  stp: number;
  /** φ at the returned step. */
  phi: number;
  /** φ' at the returned step. */
  phiPrime: number;
  /** Number of φ/φ' evaluations consumed. */
  nfev: number;
  /** True iff status indicates usable progress (convergence or warning). */
  ok: boolean;
}

/* ================================================================== */
/*  Internal state                                                     */
/* ================================================================== */

/**
 * Persistent state of the `dcsrch` state machine. Tracks the uncertainty
 * interval `[stx, sty]`, stage (ψ vs φ), and interval-shrinkage history.
 */
interface DcsrchState {
  brackt: boolean;
  stage: 1 | 2;
  ginit: number; // φ'(0)
  gtest: number; // ftol · ginit (< 0 on descent)
  finit: number; // φ(0)
  // Best (stx) and other (sty) endpoints: store φ-values.
  stx: number;
  fx: number;
  gx: number;
  sty: number;
  fy: number;
  gy: number;
  // Current admissible interval for the next trial.
  stmin: number;
  stmax: number;
  // Interval width history for the bisection safeguard.
  width: number;
  width1: number;
}

const XTRAPL = 1.1;
const XTRAPU = 4.0;
const P5 = 0.5;
const P66 = 0.66;

/**
 * Cap on consecutive non-finite φ/φ' evaluations within a single
 * NaN-bisection episode. Mirrors the Fortran v3.0 `iback ≥ 20`
 * threshold inside one line-search call (see `lnsrlb.f`); past this
 * count the line search returns `'error'` instead of bisecting the
 * interval to zero, which would otherwise be reported as silent
 * convergence at `state.stx`.
 */
const NAN_CAP = 20;

/* ================================================================== */
/*  dcstep — safeguarded cubic/quadratic step                          */
/* ================================================================== */

/**
 * One step of the Moré–Thuente interval update. Mutates `s` in place:
 *   - updates `stx/fx/gx, sty/fy/gy, brackt` per the 4-case logic,
 *   - returns the next trial step (not yet clamped to caller's cap).
 *
 * All inputs (`stp, fp, gp`) must be finite.
 */
function dcstep(
  s: DcsrchState,
  stp: number,
  fp: number,
  gp: number,
  stpmin: number,
  stpmax: number,
): number {
  const stx = s.stx;
  const fx = s.fx;
  const dx = s.gx;
  const sty = s.sty;
  const fy = s.fy;
  const dy = s.gy;
  const brackt = s.brackt;

  const sgnd = gp * Math.sign(dx);
  let stpf: number;
  let newBrackt = brackt;

  if (fp > fx) {
    // Case 1: function value rose → minimum is between stx and stp.
    const theta = 3.0 * (fx - fp) / (stp - stx) + dx + gp;
    const sScale = Math.max(Math.abs(theta), Math.abs(dx), Math.abs(gp));
    const ts = theta / sScale;
    let gamma = sScale * Math.sqrt(ts * ts - (dx / sScale) * (gp / sScale));
    if (stp < stx) gamma = -gamma;
    const p = (gamma - dx) + theta;
    const q = ((gamma - dx) + gamma) + gp;
    const r = p / q;
    const stpc = stx + r * (stp - stx);
    const stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx);
    if (Math.abs(stpc - stx) < Math.abs(stpq - stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / 2.0;
    newBrackt = true;
  } else if (sgnd < 0.0) {
    // Case 2: slope changed sign → minimum bracketed.
    const theta = 3.0 * (fx - fp) / (stp - stx) + dx + gp;
    const sScale = Math.max(Math.abs(theta), Math.abs(dx), Math.abs(gp));
    const ts = theta / sScale;
    let gamma = sScale * Math.sqrt(ts * ts - (dx / sScale) * (gp / sScale));
    if (stp > stx) gamma = -gamma;
    const p = (gamma - gp) + theta;
    const q = ((gamma - gp) + gamma) + dx;
    const r = p / q;
    const stpc = stp + r * (stx - stp);
    const stpq = stp + (gp / (gp - dx)) * (stx - stp);
    if (Math.abs(stpc - stp) > Math.abs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    newBrackt = true;
  } else if (Math.abs(gp) < Math.abs(dx)) {
    // Case 3: same sign, slope decreased in magnitude.
    const theta = 3.0 * (fx - fp) / (stp - stx) + dx + gp;
    const sScale = Math.max(Math.abs(theta), Math.abs(dx), Math.abs(gp));
    const ts = theta / sScale;
    // Cubic may have no real min → guard.
    let gamma = sScale * Math.sqrt(Math.max(0.0, ts * ts - (dx / sScale) * (gp / sScale)));
    if (stp > stx) gamma = -gamma;
    const p = (gamma - gp) + theta;
    const q = (gamma + (dx - gp)) + gamma;
    const r = p / q;
    let stpc: number;
    if (r < 0.0 && gamma !== 0.0)
      stpc = stp + r * (stx - stp);
    else if (stp > stx)
      stpc = stpmax;
    else
      stpc = stpmin;
    const stpq = stp + (gp / (gp - dx)) * (stx - stp);
    if (brackt) {
      if (Math.abs(stpc - stp) < Math.abs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      if (stp > stx)
        stpf = Math.min(stp + P66 * (sty - stp), stpf);
      else
        stpf = Math.max(stp + P66 * (sty - stp), stpf);
    } else {
      if (Math.abs(stpc - stp) > Math.abs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      stpf = Math.min(stpmax, stpf);
      stpf = Math.max(stpmin, stpf);
    }
  } else {
    // Case 4: same sign, slope did not decrease.
    if (brackt) {
      const theta = 3.0 * (fp - fy) / (sty - stp) + dy + gp;
      const sScale = Math.max(Math.abs(theta), Math.abs(dy), Math.abs(gp));
      const ts = theta / sScale;
      let gamma = sScale * Math.sqrt(ts * ts - (dy / sScale) * (gp / sScale));
      if (stp > sty) gamma = -gamma;
      const p = (gamma - gp) + theta;
      const q = ((gamma - gp) + gamma) + dy;
      const r = p / q;
      stpf = stp + r * (sty - stp);
    } else if (stp > stx)
      stpf = stpmax;
    else
      stpf = stpmin;
  }

  // Update interval endpoints. Order matters: case 2 reshuffles first.
  if (fp > fx) {
    s.sty = stp; s.fy = fp; s.gy = gp;
  } else {
    if (sgnd < 0.0) {
      s.sty = stx; s.fy = fx; s.gy = dx;
    }
    s.stx = stp; s.fx = fp; s.gx = gp;
  }
  s.brackt = newBrackt;

  return stpf;
}

/* ================================================================== */
/*  dcsrch — main line-search state machine                            */
/* ================================================================== */

/**
 * One tick of the `dcsrch` state machine. The caller alternates between
 * feeding a freshly-evaluated `(phi, phiPrime)` at `stp` and reading the
 * next trial step from the return value.
 *
 * `task === 'fg'` on the first call initialises the state; subsequently
 * the machine either returns another `'fg'` request or a terminal state.
 */
interface DcsrchTick {
  stp: number;
  task: LineSearchStatus | 'fg';
}

function dcsrch(
  s: DcsrchState,
  stp: number,
  phi: number,
  phiPrime: number,
  isFirstCall: boolean,
  ftol: number,
  gtol: number,
  xtol: number,
  stpmin: number,
  stpmax: number,
): DcsrchTick {
  // Initial call: stash descent info, set up [0, stp + XTRAPU·stp].
  if (isFirstCall) {
    s.brackt = false;
    s.stage = 1;
    s.finit = phi;
    s.ginit = phiPrime;
    s.gtest = ftol * phiPrime;
    s.width = stpmax - stpmin;
    s.width1 = 2 * s.width;
    s.stx = 0; s.fx = phi; s.gx = phiPrime;
    s.sty = 0; s.fy = phi; s.gy = phiPrime;
    s.stmin = 0;
    s.stmax = stp + XTRAPU * stp;
    return {stp, task: 'fg'};
  }

  const ftest = s.finit + stp * s.gtest;

  // Stage transition: once we have a point with sufficient decrease
  // and non-negative slope under ψ, ψ behaves like φ on the remaining
  // interval, so switch to φ-based updates.
  if (s.stage === 1 && phi <= ftest && phiPrime >= 0)
    s.stage = 2;

  // --- Termination tests, in priority order ---
  if (s.brackt && (stp <= s.stmin || stp >= s.stmax))
    return {stp: s.stx, task: 'warning_rounding'};
  if (s.brackt && (s.stmax - s.stmin) <= xtol * s.stmax)
    return {stp: s.stx, task: 'warning_xtol'};
  if (stp === stpmax && phi <= ftest && phiPrime <= s.gtest)
    return {stp, task: 'warning_stpmax'};
  if (stp === stpmin && (phi > ftest || phiPrime >= s.gtest))
    return {stp, task: 'warning_stpmin'};
  if (phi <= ftest && Math.abs(phiPrime) <= gtol * (-s.ginit))
    return {stp, task: 'converged'};

  // --- dcstep call, using ψ in stage 1 when useful ---
  let stpf: number;
  if (s.stage === 1 && phi <= s.fx && phi > ftest) {
    // Operate on ψ by subtracting the linear term; unpack φ afterwards.
    const fm = phi - stp * s.gtest;
    const fxm = s.fx - s.stx * s.gtest;
    const fym = s.fy - s.sty * s.gtest;
    const gm = phiPrime - s.gtest;
    const gxm = s.gx - s.gtest;
    const gym = s.gy - s.gtest;

    // dcstep mutates (fx/fy/gx/gy/stx/sty/brackt) in-place; run on the
    // ψ-view then translate back.
    const snap = {...s, fx: fxm, gx: gxm, fy: fym, gy: gym};
    stpf = dcstep(snap, stp, fm, gm, s.stmin, s.stmax);
    s.stx = snap.stx;
    s.sty = snap.sty;
    s.brackt = snap.brackt;
    s.fx = snap.fx + s.stx * s.gtest;
    s.fy = snap.fy + s.sty * s.gtest;
    s.gx = snap.gx + s.gtest;
    s.gy = snap.gy + s.gtest;
  } else
    stpf = dcstep(s, stp, phi, phiPrime, s.stmin, s.stmax);


  // --- Bisection safeguard: force interval to shrink by at least 2/3. ---
  if (s.brackt) {
    if (Math.abs(s.sty - s.stx) >= P66 * s.width1)
      stpf = s.stx + P5 * (s.sty - s.stx);
    s.width1 = s.width;
    s.width = Math.abs(s.sty - s.stx);
  }

  // --- Update admissible interval [stmin, stmax] for the next trial. ---
  if (s.brackt) {
    s.stmin = Math.min(s.stx, s.sty);
    s.stmax = Math.max(s.stx, s.sty);
  } else {
    s.stmin = stpf + XTRAPL * (stpf - s.stx);
    s.stmax = stpf + XTRAPU * (stpf - s.stx);
  }

  // --- Final clamp + fallback on degenerate intervals. ---
  stpf = Math.max(stpmin, Math.min(stpmax, stpf));
  if (s.brackt && (stpf <= s.stmin || stpf >= s.stmax ||
      (s.stmax - s.stmin) <= xtol * s.stmax))
    stpf = s.stx;

  return {stp: stpf, task: 'fg'};
}

/* ================================================================== */
/*  Public drivers                                                     */
/* ================================================================== */

/**
 * Run the Moré–Thuente line search starting at `alpha0` given the
 * values φ(0), φ'(0). `evalFn(α)` must return `{phi, phiPrime}` at
 * `α > 0`.
 *
 * Pre-condition: `phiPrime0 < 0` (descent direction). Returns with
 * `status === 'error'` if violated.
 */
export function runLineSearch(
  evalFn: PhiEval,
  phi0: number,
  phiPrime0: number,
  alpha0: number,
  params: LineSearchParams,
): LineSearchResult {
  const stpmin = params.stpMin ?? 1e-20;
  const stpmax = params.stpMax;
  if (!(phiPrime0 < 0) || !(alpha0 > 0) || !(stpmax >= alpha0) || !(stpmax > stpmin))
    return {status: 'error', stp: alpha0, phi: phi0, phiPrime: phiPrime0, nfev: 0, ok: false};

  const state = makeState();
  let stp = alpha0;
  let phi = phi0;
  let phiPrime = phiPrime0;
  let isFirst = true;
  let nfev = 0;

  for (let k = 0; k <= params.maxSteps; k++) {
    const tick = dcsrch(
      state, stp, phi, phiPrime, isFirst,
      params.ftol, params.gtol, params.xtol, stpmin, stpmax,
    );
    isFirst = false;
    stp = tick.stp;

    if (tick.task !== 'fg')
      return finalise(tick.task, stp, phi, phiPrime, nfev, evalFn);

    if (nfev >= params.maxSteps)
      return finalise('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);

    // Evaluate φ/φ' at the requested step, bisecting toward stx on
    // NaN/Inf trials until a finite evaluation is obtained. This keeps
    // +Infinity out of dcstep arithmetic (the paper spec §9 glosses
    // over the cancellation that +Inf causes in the cubic formula).
    // `nanCount` accumulates within a single bisection episode and
    // resets once a finite value is obtained — matches the Fortran
    // `iback` semantics scoped to one line-search call.
    let ev = evalFn(stp);
    nfev++;
    let nanCount = 0;
    while (!Number.isFinite(ev.phi) || !Number.isFinite(ev.phiPrime)) {
      nanCount++;
      if (nanCount >= NAN_CAP)
        return finalise('error', state.stx, state.fx, state.gx, nfev, evalFn);
      if (nfev >= params.maxSteps)
        return finalise('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);
      const shrunk = state.stx + 0.5 * (stp - state.stx);
      if (shrunk === stp || shrunk <= stpmin)
        return finalise('warning_rounding', state.stx, state.fx, state.gx, nfev, evalFn);
      stp = shrunk;
      ev = evalFn(stp);
      nfev++;
    }
    phi = ev.phi;
    phiPrime = ev.phiPrime;
  }
  return finalise('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);
}

/** Async variant of {@link runLineSearch}. */
export async function runLineSearchAsync(
  evalFn: PhiEvalAsync,
  phi0: number,
  phiPrime0: number,
  alpha0: number,
  params: LineSearchParams,
): Promise<LineSearchResult> {
  const stpmin = params.stpMin ?? 1e-20;
  const stpmax = params.stpMax;
  if (!(phiPrime0 < 0) || !(alpha0 > 0) || !(stpmax >= alpha0) || !(stpmax > stpmin))
    return {status: 'error', stp: alpha0, phi: phi0, phiPrime: phiPrime0, nfev: 0, ok: false};

  const state = makeState();
  let stp = alpha0;
  let phi = phi0;
  let phiPrime = phiPrime0;
  let isFirst = true;
  let nfev = 0;

  for (let k = 0; k <= params.maxSteps; k++) {
    const tick = dcsrch(
      state, stp, phi, phiPrime, isFirst,
      params.ftol, params.gtol, params.xtol, stpmin, stpmax,
    );
    isFirst = false;
    stp = tick.stp;

    if (tick.task !== 'fg')
      return finaliseAsync(tick.task, stp, phi, phiPrime, nfev, evalFn);

    if (nfev >= params.maxSteps)
      return finaliseAsync('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);

    // Same NaN/Inf bisection protocol as the sync path; see comment
    // above for `nanCount` semantics.
    let ev = await evalFn(stp);
    nfev++;
    let nanCount = 0;
    while (!Number.isFinite(ev.phi) || !Number.isFinite(ev.phiPrime)) {
      nanCount++;
      if (nanCount >= NAN_CAP)
        return finaliseAsync('error', state.stx, state.fx, state.gx, nfev, evalFn);
      if (nfev >= params.maxSteps)
        return finaliseAsync('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);
      const shrunk = state.stx + 0.5 * (stp - state.stx);
      if (shrunk === stp || shrunk <= stpmin)
        return finaliseAsync('warning_rounding', state.stx, state.fx, state.gx, nfev, evalFn);
      stp = shrunk;
      ev = await evalFn(stp);
      nfev++;
    }
    phi = ev.phi;
    phiPrime = ev.phiPrime;
  }
  return finaliseAsync('max_steps_reached', stp, phi, phiPrime, nfev, evalFn);
}

/* ================================================================== */
/*  Helpers                                                            */
/* ================================================================== */

function makeState(): DcsrchState {
  return {
    brackt: false, stage: 1,
    ginit: 0, gtest: 0, finit: 0,
    stx: 0, fx: 0, gx: 0,
    sty: 0, fy: 0, gy: 0,
    stmin: 0, stmax: 0,
    width: 0, width1: 0,
  };
}

function finalise(
  status: LineSearchStatus,
  stp: number,
  phi: number,
  phiPrime: number,
  nfev: number,
  evalFn: PhiEval,
): LineSearchResult {
  // If dcsrch returned a fallback stp (e.g. state.stx on warning) the
  // stored phi/phiPrime may not correspond to that stp — re-evaluate.
  if (status !== 'converged' && status !== 'max_steps_reached') {
    const ev = evalFn(stp);
    if (Number.isFinite(ev.phi) && Number.isFinite(ev.phiPrime)) {
      phi = ev.phi;
      phiPrime = ev.phiPrime;
      nfev++;
    }
  }
  const ok = status !== 'error' && status !== 'max_steps_reached';
  return {status, stp, phi, phiPrime, nfev, ok};
}

async function finaliseAsync(
  status: LineSearchStatus,
  stp: number,
  phi: number,
  phiPrime: number,
  nfev: number,
  evalFn: PhiEvalAsync,
): Promise<LineSearchResult> {
  if (status !== 'converged' && status !== 'max_steps_reached') {
    const ev = await evalFn(stp);
    if (Number.isFinite(ev.phi) && Number.isFinite(ev.phiPrime)) {
      phi = ev.phi;
      phiPrime = ev.phiPrime;
      nfev++;
    }
  }
  const ok = status !== 'error' && status !== 'max_steps_reached';
  return {status, stp, phi, phiPrime, nfev, ok};
}
