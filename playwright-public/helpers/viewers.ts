/**
 * Playwright helpers shared by spec files in the public CI mirror.
 *
 * The dev `TestTrack/helpers/viewers.ts` carries the full viewer/legend helper
 * surface; the public mirror only needs `finishSpec`, so this is the trimmed
 * port. Keep the same export name/signature so the mirrored specs import it
 * unchanged.
 */

import {Page} from '@playwright/test';
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

/**
 * Canonical Cleanup step: `grok.shell.closeAll()` + a 500ms settle. Optionally
 * clears the categorical color-coding tags on the active table's `Stereo
 * Category` column (used by chem-flavoured charts specs). Mirrors the dev
 * `TestTrack/helpers/viewers.ts` helper so the mirrored specs import it
 * unchanged.
 */
export async function cleanupShell(
  page: Page,
  opts: {clearStereoCategoryColorCoding?: boolean} = {},
): Promise<void> {
  await page.evaluate((clearColors) => {
    if (clearColors) {
      try {
        const col = (window as any).grok.shell.tv?.dataFrame.col('Stereo Category');
        if (col) {
          delete col.tags['.color-coding-categorical'];
          delete col.tags['.color-coding-type'];
        }
      } catch (_) {}
    }
    (window as any).grok.shell.closeAll();
  }, opts.clearStereoCategoryColorCoding ?? false);
  await page.waitForTimeout(500);
}
