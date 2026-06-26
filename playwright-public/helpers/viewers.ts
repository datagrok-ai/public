/**
 * Playwright helpers shared by spec files in the public CI mirror.
 *
 * The dev `TestTrack/helpers/viewers.ts` carries the full viewer/legend helper
 * surface; the public mirror only needs `finishSpec`, so this is the trimmed
 * port. Keep the same export name/signature so the mirrored specs import it
 * unchanged.
 */

import {stepErrors, StepError} from '../spec-login';

/**
 * Trailing soft-step assertion (throws if any softStep failed). Reads from the
 * shared `stepErrors` array exported by spec-login. Pass a non-default `prefix`
 * only if a spec wants a different message header.
 */
export function finishSpec(prefix = 'Step failures'): void {
  if (stepErrors.length === 0) return;
  const summary = stepErrors.map((e: StepError) => `- ${e.step}: ${e.error}`).join('\n');
  throw new Error(`${prefix}:\n${summary}`);
}
