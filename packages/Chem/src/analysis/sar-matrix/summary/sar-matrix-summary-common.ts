import {RoleSummary} from './sar-matrix-role-fit';
import {SwapPool, TRUST_R2} from './sar-matrix-summary-data';

/** Component rows shown before the rest fold away, so the swaps stay on the first screen. */
export const FINDING_ROWS = 4;

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

export interface Answer {
  /** Absent where the drawing is the answer. */
  answer?: string;
  negative: string;
  /** Drawn, since truncated SMILES of related analogs read alike. */
  art?: HTMLElement;
}

export function shortSmiles(smiles: string): string {
  return smiles.length > 20 ? `${smiles.slice(0, 19)}…` : smiles;
}

export function roleRanks(role: RoleSummary): boolean {
  return role.levels.length >= 2 && (role.repeat === null || role.repeat >= TRUST_R2);
}
