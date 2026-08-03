import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// playwright-public hosts the core/platform E2E suites. All general config lives in
// the shared base (@datagrok-libraries/test/src/playwright/base-config); here we only
// set what is specific to this run dir.
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
  // Kill-switch so a degraded run stops itself and still writes its JSON report + CSV
  // instead of being SIGKILL-aborted at the 240-min Jenkins stage timeout (which loses
  // all results). A healthy no-retry run is ~65 min; 150 min leaves >2x headroom for a
  // slow/contended stand while capping waste well under the stage kill.
  globalTimeout: 150 * 60 * 1000,
});
