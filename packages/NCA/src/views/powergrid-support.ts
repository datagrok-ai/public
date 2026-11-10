/**
 * Does the installed PowerGrid draw the `Violated` reason pills readably?
 *
 * The pills are drawn by PowerGrid's `Tags` cell renderer. Released PowerGrid ≤ 1.8.4 fills each
 * pill without `beginPath()`, so the grid's cell rectangle rides along in the path and the whole
 * cell takes the last pill's colour — only the last reason stays visible (measured 2026-10-19). The
 * fix (GROK-20888) sits in PowerGrid's unreleased `v.next`, so every release AFTER 1.8.4 carries it.
 *
 * Checked once, when the app starts ({@link detectPowerGridSupport} from `ncaApp` /
 * `ncaAppTreeBrowser`), and cached; any later reader gets the same answer. The answer only sets the
 * DEFAULT of the Settings switch — the scientist can still turn the column on or off.
 */

import * as DG from 'datagrok-api/dg';

/** The first PowerGrid release that can carry the GROK-20888 fix. */
export const PILLS_MIN_POWERGRID = '1.8.5';

/** A debug build (`<version>.X-<hash>`) of this version or later is accepted too: it is built from
 *  the current source, which carries the fix (user-directed 2026-10-19). */
export const PILLS_MIN_POWERGRID_DEBUG = '1.8.4';

/** True for a debug build: `1.8.4.X-44ea6cd6`. */
function isDebugBuild(version: string): boolean {
  return /\.X(-|$)/.test(version);
}

/** Does this PowerGrid version draw the pills readably? A release ≥ {@link PILLS_MIN_POWERGRID}, or a
 *  debug build of ≥ {@link PILLS_MIN_POWERGRID_DEBUG}. */
export function powerGridDrawsPills(version: string): boolean {
  return versionAtLeast(version, PILLS_MIN_POWERGRID) ||
    (isDebugBuild(version) && versionAtLeast(version, PILLS_MIN_POWERGRID_DEBUG));
}

export interface PowerGridSupport {
  /** The PowerGrid version this session loaded, or null when PowerGrid is not installed. */
  readonly version: string | null;
  /** True when that version draws the pills readably — the default of the `Violated` switch. */
  readonly pills: boolean;
}

/** The leading numeric parts of a version: `1.8.4.X-44ea6cd6` → [1, 8, 4]. */
function numericParts(version: string): number[] {
  const out: number[] = [];
  for (const part of version.split('.')) {
    if (!/^\d+$/.test(part)) break;
    out.push(Number(part));
  }
  return out;
}

/** True when `version` ≥ `min`, comparing the numeric parts only. */
export function versionAtLeast(version: string, min: string): boolean {
  const a = numericParts(version);
  const b = numericParts(min);
  if (a.length === 0) return false;
  for (let i = 0; i < Math.max(a.length, b.length); i++) {
    const x = a[i] ?? 0;
    const y = b[i] ?? 0;
    if (x !== y) return x > y;
  }
  return true;
}

/** The version of the PowerGrid package that provides the `Tags` renderer in this session. */
function loadedPowerGridVersion(): string | null {
  const fromPackage = DG.Func.find({package: 'PowerGrid'});
  const tags = fromPackage.find((f) => f.name === 'Tags') ?? fromPackage[0];
  const version = tags?.package?.version;
  return typeof version === 'string' && version.length > 0 ? version : null;
}

let cached: PowerGridSupport | null = null;

/** The PowerGrid check — run at app start, cached for the session. Never throws: an API failure
 *  reads as "not installed", which leaves the column off (the safe side). */
export function detectPowerGridSupport(): PowerGridSupport {
  if (cached !== null) return cached;
  let version: string | null = null;
  try {
    version = loadedPowerGridVersion();
  } catch (e) {
    console.warn('[NCA] PowerGrid version check failed — the Violated column starts hidden.', e);
  }
  cached = {version, pills: version !== null && powerGridDrawsPills(version)};
  return cached;
}
