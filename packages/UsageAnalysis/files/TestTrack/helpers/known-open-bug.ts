/**
 * known-open-bug — self-flipping wrapper for a STILL-OPEN bug's regression
 * assertion.
 *
 * A scenario that guards an open bug keeps the REAL, strict assertion on the
 * desired (currently-broken) behaviour, wrapped so an expected reproduction is
 * green and an unexpected fix is loud:
 *
 *   - bug open + assertion FAILS (bug reproduces) → green; logs
 *     `[KNOWN_BUG_REPRODUCED:<id>]` (Gate B folds this into the
 *     KNOWN_BUG_REPRODUCED verdict — no repair loop).
 *   - bug open + assertion PASSES (bug looks fixed) → loud throw
 *     `[KNOWN_BUG_FIXED:<id>]` (Gate B failure key B-KNOWN-BUG-FIXED-RESTORE):
 *     the operator sets `related_bugs[].status: fixed` / `fixed_in` and replaces
 *     this wrapper with a plain hard `expect`.
 *
 * Strictly better than a `console.warn` downgrade: the real `expect` still runs
 * (audit-visible, not a weakening) and the convention self-flips when devs fix
 * the bug. Use ONLY for an explicit `related_bugs[].status: open` entry; a
 * fixed / bare-string bug ref uses a plain hard `expect`.
 */
export async function knownOpenBug(
  bugId: string,
  assertion: () => void | Promise<void>,
): Promise<void> {
  try {
    await assertion();
  } catch {
    console.log(`[KNOWN_BUG_REPRODUCED:${bugId}]`);
    return;
  }
  throw new Error(
    `[KNOWN_BUG_FIXED:${bugId}] assertion now passes — set related_bugs[].status: ` +
    `fixed + fixed_in and replace knownOpenBug with a plain hard expect`,
  );
}
