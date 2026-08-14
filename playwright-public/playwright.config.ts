import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// playwright-public hosts the core/platform E2E suites. All general config lives in
// the shared base (@datagrok-libraries/test/src/playwright/base-config); here we only
// set what is specific to this run dir.
// Refuse to run without an explicit target. baseConfig falls back to
// https://dev.datagrok.ai when DATAGROK_URL is unset, so a harness that forgot to export
// it would point the whole suite at a live environment — creating projects, editing
// shares and deleting entities there — and nothing in the output would say so. Tests
// belong on the ephemeral CI stand; failing loudly is the only safe default.
if (!process.env.DATAGROK_URL)
  throw new Error('DATAGROK_URL is not set. Run via `grok test --host <ci-stand>`; ' +
    'this suite must never fall back to a deployed environment.');

export default defineConfig({
  ...baseConfig,
  testDir: '.',
  // `helpers/` still holds the headed-only `session-helpers.test.ts` (manual second
  // user, not CI-runnable). Exclude the folder from discovery.
  testIgnore: ['**/helpers/**'],
  // This is a nightly *measurement* run (feeds dashboards, not a merge gate) of ~360
  // sequential specs. The CI pipeline passes `grok test --no-retry`, but that flag is
  // silently dropped before it reaches Playwright (minimist parses `--no-retry` as
  // `{retry:false}`, so the runner's `args['no-retry']` check is never true), leaving
  // baseConfig's `retries: CI ? 1 : 0` in force. On the shared Build-Deploy stand many
  // specs fail for environmental reasons (stand pollution + CPU contention), and each
  // failure was being re-run once — ~30% wasted wall time on a good run and enough to
  // push a degraded run past the 240-min stage kill (SIGKILL → no report at all).
  // Retries recovered 2 of 361 specs, so force them off here for this suite.
  retries: 0,
  // ~360 specs used to run strictly serially (~150 min) — the whole reason a Build-Deploy
  // Playwright stage sat for three hours. `fullyParallel` stays false, so ordering WITHIN
  // a file is untouched and only separate files run side by side.
  //
  // Measured on the 32-vCPU agent: 4 workers = 61 min / 55 failures, 8 = 47 min / 59.
  // The extra four land in browse, connections and projects, and projects tracks
  // concurrency directly (2 serial, 5 at four, 6 at eight), so those files share
  // server-side state. Fourteen minutes is not worth four phantom failures in a suite
  // whose whole job is signal; the way to 8 is isolating that state, not this knob.
  // PLAYWRIGHT_WORKERS still overrides per run.
  workers: Number(process.env.PLAYWRIGHT_WORKERS ?? 4),
  // Kill-switch so a degraded run stops itself and still writes its JSON report + CSV
  // instead of being SIGKILL-aborted at the 240-min Jenkins stage timeout (which loses
  // all results). A healthy no-retry run is ~65 min; 150 min leaves >2x headroom for a
  // slow/contended stand while capping waste well under the stage kill.
  globalTimeout: 150 * 60 * 1000,
});
