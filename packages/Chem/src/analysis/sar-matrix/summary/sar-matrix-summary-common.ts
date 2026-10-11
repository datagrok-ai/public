/* What every Summary section shares: the section names, the thresholds the rendering reads and the
   small pure helpers that format and orient what the collector produced. */

import {RoleSummary} from './sar-matrix-role-fit';
import {SwapPool, TRUST_R2} from './sar-matrix-summary-data';

/** Component rows shown before the rest fold away, so the swaps stay on the first screen. */
export const FINDING_ROWS = 4;

/** Spread difference below which two components tie: under the two decimals printed. */
export const ROLE_SPREAD_TIE = 0.01;

export const PANE_OVERVIEW = 'Overview';
export const PANE_EFFECTS = 'Effects';
export const PANE_MAKING = 'Worth making';
/** The Effects view of what was measured inside each series, as against the cross-table fit. */
export const PANE_SERIES = 'Measured in series';
export const PANES = [PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING];

/** Index-aligned with the five picks of Start here. */
export const REASON_GLYPHS = ['n', '↔', '★', '✚', '✓'];
export const REASON_WORDS = ['most measured compounds', 'widest measured range',
  'holds the best measured compound', 'most predictions worth making', 'best-validated fit'];

/** A pool read in the direction that improves potency; the reverse is the same data with signs flipped. */
export function orient(pool: SwapPool): {forward: boolean, from: string, to: string, worst: number,
  widest: number, mean: number, up: number} {
  const forward = pool.min >= -pool.max;
  return {
    forward,
    from: forward ? pool.from : pool.to,
    to: forward ? pool.to : pool.from,
    worst: forward ? pool.min : -pool.max,
    widest: forward ? pool.max : -pool.min,
    mean: (forward ? pool.sum : -pool.sum) / pool.n,
    up: forward ? pool.nUp : pool.n - pool.nUp,
  };
}

/** One answer tile: what it is about, the conclusion, and the clause that limits it. */
export interface Answer {
  answer: string;
  negative: string;
  /** Drawn, since truncated SMILES of related analogs read alike. */
  art?: HTMLElement;
}

/** A fragment SMILES cut to fit a one-line descriptor; the full string stays in the tooltip. */
export function shortSmiles(smiles: string): string {
  return smiles.length > 20 ? `${smiles.slice(0, 19)}…` : smiles;
}

/** The anchor of a component's swap card on Effects. */
export function swapAnchor(role: string): string {
  return `swaps-${role}`;
}

/** A blank value is the unsubstituted parent, a level like any other, so it is named. */
export function roleValueName(value: string): string {
  return value === '' ? 'nothing at this position' : value;
}

/** Whether a component's values may be ranked: two readable ones, and offsets that repeat across halves. */
export function roleRanks(role: RoleSummary): boolean {
  return role.levels.length >= 2 && (role.repeat === null || role.repeat >= TRUST_R2);
}
