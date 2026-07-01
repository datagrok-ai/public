/* eslint-disable */
// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./line-search.ts
import {LineSearchParams, LineSearchResult, LineSearchStatus, NAN_CAP, dcsrch, makeState} from './line-search';

export function runLineSearch(
  evalFn: (alpha: number) => {phi: number; phiPrime: number},
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

function finalise(
  status: LineSearchStatus,
  stp: number,
  phi: number,
  phiPrime: number,
  nfev: number,
  evalFn: (alpha: number) => {phi: number; phiPrime: number},
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
