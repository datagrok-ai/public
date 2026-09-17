/* How long an assertion is given to come true. Every check in the library goes through this
   `expect`, so the budget can be narrowed for a stretch of a run: a `@known-failure` scenario
   states a defect the product still has, and letting each of its assertions spend the suite's
   15 s expect budget proving what the tag already says costs a minute of a run doing nothing. */
import {expect as playwright} from '@playwright/test';

/** What an assertion inside a `@known-failure` scenario waits: long enough for a fixed bug to pass
 * on the second poll — the other half of what the tag checks — and short enough that a run is not
 * spent on it. */
export const KNOWN_FAILURE_MS = 3000;

let current = playwright;
let narrowed = false;

export const expect: typeof playwright = new Proxy(playwright, {
  apply: (_t, _this, args: unknown[]) => (current as (...a: unknown[]) => unknown)(...args),
  get: (_t, prop) => Reflect.get(current, prop),
}) as typeof playwright;

/** Runs a scenario whose failure is expected, with the narrowed budget. */
export async function whileExpectedToFail<T>(body: () => Promise<T>): Promise<T> {
  current = playwright.configure({timeout: KNOWN_FAILURE_MS});
  narrowed = true;
  try {
    return await body();
  }
  finally {
    current = playwright;
    narrowed = false;
  }
}

/** The budget a check names for itself, narrowed the same way — `expect.configure` only reaches
 * the checks that name none. */
export function pollMs(ms: number): number {
  return narrowed ? Math.min(ms, KNOWN_FAILURE_MS) : ms;
}
