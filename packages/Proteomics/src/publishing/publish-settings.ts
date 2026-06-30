/**
 * Configurable destination for shared review snapshots.
 *
 * Both values are package settings (declared in `package.json` under `properties`, edited
 * from the package's **Settings** panel in the platform). They are read at share time, so
 * an admin can retarget where reviews land without a code change.
 *
 * What is and isn't configurable, and why:
 *   - The **umbrella Space** ({@link reviewSpaceName}) and the **name prefix**
 *     ({@link reviewNamePrefix}) are free text.
 *   - The structural suffix `-<target-slug>-v<version>-<date>` is NOT configurable: the
 *     republish flow finds the prior version by matching `"<prefix>-<slug>-v"` and parsing
 *     the `-v<n>-` segment (see `findPriorShare` in `publish-state.ts`). Letting the suffix
 *     vary would break version detection. So we expose the prefix, keep the skeleton fixed.
 *
 * Note: changing the prefix after shares already exist orphans the old-prefix versions from
 * republish detection (a fresh share starts at v1 again). That's the expected cost of
 * retargeting; the defaults reproduce the original hard-coded behaviour.
 */
import {_package} from '../package';

export const DEFAULT_REVIEW_SPACE = 'Proteomics-Reviews';
export const DEFAULT_REVIEW_NAME_PREFIX = 'Proteomics-Review';

/** Reads a string package setting, falling back when unset/blank/not-yet-loaded. */
function readSetting(key: string, fallback: string): string {
  try {
    const v = (_package?.settings as Record<string, unknown> | undefined)?.[key];
    if (typeof v === 'string' && v.trim().length > 0)
      return v.trim();
  } catch { /* settings not ready — use fallback */ }
  return fallback;
}

/** Umbrella Space the per-target review Spaces live under. Setting: `reviewSpaceName`. */
export function reviewSpaceName(): string {
  return readSetting('reviewSpaceName', DEFAULT_REVIEW_SPACE);
}

/** Prefix for the per-target child Space and the published project name.
 * Setting: `reviewNamePrefix`. The `-<slug>-v<version>-<date>` suffix stays fixed. */
export function reviewNamePrefix(): string {
  return readSetting('reviewNamePrefix', DEFAULT_REVIEW_NAME_PREFIX);
}

/**
 * Whether to re-open each published dashboard and assert it survives a reload
 * (the heavy round-trip check). Default ON — keep it for client deliverables;
 * turn it off via the `verifyPublishedDashboard` package setting for faster demo
 * shares. The reviewer-access verify-and-rollback gate is separate and always runs.
 */
export function verifyPublishedDashboard(): boolean {
  try {
    const v = (_package?.settings as Record<string, unknown> | undefined)?.['verifyPublishedDashboard'];
    if (v === false || v === 'false') return false;
    if (v === true || v === 'true') return true;
  } catch { /* settings not ready — use default */ }
  return true;
}
