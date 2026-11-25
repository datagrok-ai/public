/**
 * Route-code contract (sci-comp-ip-route, peer review R2 N-1): every `RouteCode`
 * 0…6 is computed; a value outside the set (7, NaN, a non-integer, negative —
 * reachable from untyped callers) is REJECTED by `computeNca`, `aucInterval` and
 * `sparseAuc`, never computed as extravascular. Before 0.13.0 such a value fell
 * through the IV checks and silently took the extravascular `(0, 0)` path.
 */

import {computeNca} from '../compute-nca';
import {aucInterval} from '../auc-interval';
import {sparseAuc} from '../sparse';
import {
  ROUTE_IV_BOLUS, ROUTE_IV_INFUSION, ROUTE_PO, ROUTE_SC, ROUTE_IM, ROUTE_OTHER, ROUTE_IP,
  isKnownRouteCode,
} from '../types';
import type {ProfileInputs, NcaRules, RouteCode} from '../types';
import {mouseInputs} from './sparse-fixture-io';

const RULES: NcaRules = {
  aucMethod: 'linear-up-log-down',
  blq: {
    preFirstMeasurable: 'set-zero',
    embedded: 'set-zero',
    afterLast: 'set-zero',
    consecutiveAfterLast: 'set-zero',
  },
  lambdaZ: {
    mode: 'auto-best-fit',
    minPoints: 3,
    minRSquared: 0.85,
    excludeCmax: true,
    adjRSquaredFactor: 1e-4,
  },
  extrapWarnPct: 20,
  extrapErrorPct: 50,
  extrapWarnPctAumc: 20,
  compensatedSummation: false,
};

const TIME = [0.25, 0.5, 1, 2, 4, 6, 8, 12];
const CONC = [1.5, 2.4, 3.0, 1.8, 0.9, 0.5, 0.3, 0.1];

function inputs(route: number): ProfileInputs {
  return {
    time: Float64Array.from(TIME),
    conc: Float64Array.from(CONC),
    blqMask: new Uint8Array(TIME.length),
    lloq: 0.01,
    dose: 2.5,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    route: route as RouteCode, // untyped caller: the cast is the point of the test
    infusionDuration: null,
    bodyWeight: null,
  };
}

const KNOWN: ReadonlyArray<RouteCode> = [
  ROUTE_IV_BOLUS, ROUTE_IV_INFUSION, ROUTE_PO, ROUTE_SC, ROUTE_IM, ROUTE_OTHER, ROUTE_IP,
];
const UNKNOWN: ReadonlyArray<number> = [7, -1, 2.5, NaN, Infinity];

describe('isKnownRouteCode', () => {
  it('accepts exactly the integers 0…6', () => {
    expect(KNOWN).toEqual([0, 1, 2, 3, 4, 5, 6]);
    for (const c of KNOWN) expect(isKnownRouteCode(c)).toBe(true);
    for (const c of UNKNOWN) expect(isKnownRouteCode(c)).toBe(false);
  });
});

describe('computeNca — route-code contract', () => {
  it.each(KNOWN)('route %d is computed (status ok, no UNKNOWN_ROUTE_CODE)', (route) => {
    const r = computeNca(inputs(route), RULES);
    expect(r.status).toBe('ok');
    expect(r.provenance.warnings.some((w) => w.code === 'UNKNOWN_ROUTE_CODE')).toBe(false);
  });

  it.each(UNKNOWN)('route %p is rejected: failed, NaN values, one error warning', (route) => {
    const r = computeNca(inputs(route), RULES);
    expect(r.status).toBe('failed');
    expect(Object.values(r.values).every((v) => Number.isNaN(v))).toBe(true);
    expect(r.provenance.warnings).toHaveLength(1);
    expect(r.provenance.warnings[0].code).toBe('UNKNOWN_ROUTE_CODE');
    expect(r.provenance.warnings[0].severity).toBe('error');
    expect(r.provenance.lambdaZ).toBeNull();
    expect(r.provenance.c0).toBeNull();
  });
});

describe('aucInterval — route-code contract', () => {
  it.each(KNOWN)('route %d gives a finite AUC', (route) => {
    expect(Number.isFinite(aucInterval(inputs(route), RULES, 0, Infinity).value)).toBe(true);
  });

  it.each(UNKNOWN)('route %p is rejected (NaN)', (route) => {
    const r = aucInterval(inputs(route), RULES, 0, Infinity);
    expect(r.value).toBeNaN();
    expect(r.tLast).toBeNaN();
  });
});

describe('sparseAuc — route-code contract', () => {
  const {input} = mouseInputs()[0];

  it.each(KNOWN)('route %d gives a finite AUC', (route) => {
    const r = sparseAuc(input, {route});
    expect(Number.isFinite(r.auc)).toBe(true);
    expect(r.warnings.some((w) => w.code === 'SPARSE_UNKNOWN_ROUTE_CODE')).toBe(false);
  });

  it.each(UNKNOWN)('route %p is rejected: NaN AUC/SE/df/CI, composite still reported', (route) => {
    const r = sparseAuc(input, {route: route as RouteCode});
    expect(r.auc).toBeNaN();
    expect(r.se).toBeNaN();
    expect(r.df).toBeNaN();
    expect(r.ci[0]).toBeNaN();
    expect(r.ci[1]).toBeNaN();
    expect(r.warnings.map((w) => w.code)).toEqual(['SPARSE_UNKNOWN_ROUTE_CODE']);
    expect(r.composite).toEqual(sparseAuc(input).composite);
    expect(r.topology).toBe(sparseAuc(input).topology);
  });
});
