/* What every Summary segment shares: the segment names, the layout breakpoints, the thresholds the
   rendering reads and the small pure helpers that format and orient what the collector produced. */

import {RoleSummary} from './sar-matrix-role-fit';
import {SwapPool, TRUST_R2} from './sar-matrix-summary-data';

export const TRUST_LIST_MAX = 10;
/** Dock widths and height below which the tab drops marks for their text, in px. */
export const NARROW_PX = 720;
export const XNARROW_PX = 560;
export const SHORT_PX = 340;
/** Component rows the landing band shows before folding the rest away, so the swap row under them and
 *  the band under that stay on the screen the reader lands on. */
export const FINDING_ROWS = 4;

/** Shortest a strip square's directional fill may draw and still read as a value rather than the
 *  mid-rule itself. */
export const STRIP_MIN_FILL = 3;
/** Multiples of a series' own leave-one-out error a tick run shows before a count stops reading at a
 *  glance. */
export const GAIN_TICKS = 4;
/** Difference in noise-corrected spread below which two component columns are one answer with two
 *  names: offsets print at two decimals, so anything under this is not on screen at all. */
export const ROLE_SPREAD_TIE = 0.01;

/** Said on the Overview and again on the core card, which a reader can arrive at either way round. */
export const CORES_NOT_COMPARABLE = 'Cores are not comparable across series here — a core is one row of one ' +
  'matrix and recurs only inside its own fold lineage. Each series\' best core is in its Start-here ' +
  'expand.';

export const PANE_OVERVIEW = 'Overview';
export const PANE_EFFECTS = 'Effects';
export const PANE_MAKING = 'Worth making';
export const PANE_METHOD = 'Method';
/** The Effects segment's last tab: what is read off measured compounds inside each series — which
 *  substituent came first where, which core scored best, which swap was actually made — as against the
 *  component tabs, which are one fitted model over the whole table. */
export const PANE_SERIES = 'Measured in series';
export const PANES = [PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING, PANE_METHOD];

/** Fixed slots of the reason lane, in this order on every row — the lane is only comparable down a
 *  column if slot 3 means the same thing everywhere. */
export const REASON_GLYPHS = ['n', '↔', '★', '✚', '✓'];
export const REASON_WORDS = ['most measured compounds', 'widest measured range',
  'holds the best measured compound', 'most predictions worth making', 'best-validated fit'];

/** A pool read in the direction that improves potency. The pool is keyed on the lexically ordered
 *  fragment pair, so the improving direction is whichever end of its range survives; the reverse swap
 *  is the same measurements with every sign flipped. */
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

/** One answer tile's content: the structures the conclusion is about, the conclusion, and the one clause
 *  that stops the conclusion being over-read. */
export interface Answer {
  answer: string;
  negative: string;
  /** Drawn above the conclusion. A substituent written out is truncated to fit the tile, and three
   *  analogs of one scaffold truncate to the same string. */
  art?: HTMLElement;
}

/** A fragment SMILES cut to fit a one-line descriptor; the full string stays in the tooltip. */
export function shortSmiles(smiles: string): string {
  return smiles.length > 20 ? `${smiles.slice(0, 19)}…` : smiles;
}

/** A blank role value is the unsubstituted parent, which is a level like any other — rendered as a
 *  name rather than as the empty string it is stored as.
 *
 *  Not shortened: a fixed character cap renders three analogs of one scaffold as the same string, since
 *  what distinguishes them is past the cap. A row ellipsizes in CSS instead, to whatever width it has. */
/** Where a component's swap card sits on the Effects segment, so the band row for that component scrolls
 *  to the card about that component rather than to whichever one was built first. */
export function swapAnchor(role: string): string {
  return `swaps-${role}`;
}

export function roleValueName(value: string): string {
  return value === '' ? 'nothing at this position' : value;
}

/** Whether this role's own values may be put in an order at all: two of them have to be readable, and
 *  a split of credit that does not survive refitting each half of the table is not this role's. */
export function roleRanks(role: RoleSummary): boolean {
  return role.levels.length >= 2 && (role.repeat === null || role.repeat >= TRUST_R2);
}
