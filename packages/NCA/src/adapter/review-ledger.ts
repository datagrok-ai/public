/**
 * Review-decision ledger operations (nca-workflow-stages F-S2.2; D11, D12, D13; 21 CFR 11.10(e)).
 *
 * One {@link ReviewLedger} per review unit holds the CURRENT decision (or `null` = `To review`)
 * and every superseded decision, oldest first. A decision is NEVER deleted: reset, auto-reset and
 * replacement all move it to `history` with the reason and time. The analysis entry persists the
 * ledgers (v15 `reviewDecisions`), so the history survives Save → reopen (AC-S2.2.1) — the
 * in-memory audit log cannot carry it (research G1).
 *
 * Pure: every function returns a NEW map; the caller (the controller) owns the state and the
 * audit-log appends.
 */

import type {ReviewDecisionRecord, ReviewLedger, SupersededBy} from './study-schema';

export type LedgerMap = ReadonlyMap<string, ReviewLedger>;

/** Build the working map from the persisted list (drops malformed / duplicate keys, first wins). */
export function ledgerMapOf(list: ReadonlyArray<ReviewLedger> | undefined): Map<string, ReviewLedger> {
  const m = new Map<string, ReviewLedger>();
  for (const l of list ?? []) if (typeof l?.unitKey === 'string' && !m.has(l.unitKey)) m.set(l.unitKey, l);
  return m;
}

/** The persisted form (only ledgers with something to keep). */
export function ledgerListOf(m: LedgerMap): ReviewLedger[] {
  return [...m.values()].filter((l) => l.current !== null || l.history.length > 0);
}

function supersede(l: ReviewLedger, by: SupersededBy, at: string): ReviewLedger {
  if (l.current === null) return l;
  return {
    unitKey: l.unitKey, current: null,
    history: [...l.history, {...l.current, supersededBy: by, supersededAt: at}],
  };
}

/** Record a new decision on a unit. Any current decision moves to history as `replaced`. */
export function recordDecision(
  m: LedgerMap, unitKey: string, rec: ReviewDecisionRecord,
): Map<string, ReviewLedger> {
  const out = new Map(m);
  const prev = out.get(unitKey) ?? {unitKey, current: null, history: []};
  out.set(unitKey, {...supersede(prev, 'replaced', rec.at), current: rec});
  return out;
}

/** Move a unit's current decision to history (`reset` by the analyst, or `replaced` by a λz
 *  point edit — D13). Returns the new map and whether anything changed. */
export function supersedeDecision(
  m: LedgerMap, unitKey: string, by: SupersededBy, at: string,
): {map: Map<string, ReviewLedger>; changed: boolean} {
  const l = m.get(unitKey);
  if (l === undefined || l.current === null) return {map: new Map(m), changed: false};
  const out = new Map(m);
  out.set(unitKey, supersede(l, by, at));
  return {map: out, changed: true};
}

/**
 * Auto-reset every current decision whose fingerprint no longer matches what is on screen
 * (D13: the fit, the source data, or the mapping changed since the analyst decided), or whose
 * unit no longer needs review / no longer exists. Returns the new map and the reset unit keys —
 * the count drives the "N review decisions reset" notice (F-S3.6).
 *
 * @param fingerprints  current decision fingerprint per unit that NEEDS REVIEW. A unit absent
 *                      from the map (now Accepted, Not applicable, or gone) cannot keep a decision.
 * @param inScope       only ledgers whose key passes are reconciled (others are left as they are).
 */
export function reconcileLedgers(
  m: LedgerMap, fingerprints: ReadonlyMap<string, string>, at: string,
  inScope: (unitKey: string) => boolean = () => true,
): {map: Map<string, ReviewLedger>; resetKeys: string[]} {
  const out = new Map(m);
  const resetKeys: string[] = [];
  for (const [key, l] of m) {
    if (l.current === null || !inScope(key)) continue;
    const fp = fingerprints.get(key);
    if (fp !== undefined && fp === l.current.fingerprint) continue;
    out.set(key, supersede(l, 'auto-reset', at));
    resetKeys.push(key);
  }
  return {map: out, resetKeys};
}
