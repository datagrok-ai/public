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
 *   - The structural suffix `-<project-slug>-v<version>-<date>` is NOT configurable: the
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

/** Umbrella Space the per-project review Spaces live under. Setting: `reviewSpaceName`. */
export function reviewSpaceName(): string {
  return readSetting('reviewSpaceName', DEFAULT_REVIEW_SPACE);
}

/** Prefix for the per-project child Space and the published project name.
 * Setting: `reviewNamePrefix`. The `-<slug>-v<version>-<date>` suffix stays fixed. */
export function reviewNamePrefix(): string {
  return readSetting('reviewNamePrefix', DEFAULT_REVIEW_NAME_PREFIX);
}

/**
 * Name of the team pre-selected in the Share for Review dialog's "Share with
 * team" picker. Setting: `defaultReviewerGroup`. Empty (the default) means "no
 * preference" — the dialog falls back to the first available team. A stale value
 * (team renamed/removed) is ignored the same way.
 */
export function defaultReviewerGroup(): string {
  return readSetting('defaultReviewerGroup', '');
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

/**
 * The controlled list of project names an analyst may share under. Maintained by
 * package administrators via the `projectVocabulary` package setting — a
 * comma-separated string they edit in the package Settings panel (add / modify /
 * remove names); the Share for Review dialog offers exactly these (no free text).
 *
 * Parses the setting robustly: a comma / newline-separated string (the string
 * property), or a JS array (future-proofing if the property is ever promoted to
 * a `string_list`). Returns a de-duplicated, trimmed, order-preserving list;
 * `[]` when unset or not-yet-loaded, which the dialog treats as "ask an admin".
 */
export function getProjectVocabulary(): string[] {
  let raw: unknown;
  try {
    raw = (_package?.settings as Record<string, unknown> | undefined)?.['projectVocabulary'];
  } catch { /* settings not ready */ return []; }
  return parseProjectVocabulary(raw);
}

/**
 * Pure parser for the `projectVocabulary` setting value — extracted so it can be
 * unit-tested without a live package. Accepts a comma / newline-separated string
 * (the string property) or a JS array (future `string_list`). Returns a
 * de-duplicated, trimmed, order-preserving list.
 */
export function parseProjectVocabulary(raw: unknown): string[] {
  let items: string[] = [];
  if (Array.isArray(raw))
    items = raw.map((v) => String(v));
  else if (typeof raw === 'string')
    items = raw.split(/[\r\n,]+/);

  const seen = new Set<string>();
  const out: string[] = [];
  for (const it of items) {
    const t = it.trim();
    if (t.length > 0 && !seen.has(t)) {
      seen.add(t);
      out.push(t);
    }
  }
  return out;
}
