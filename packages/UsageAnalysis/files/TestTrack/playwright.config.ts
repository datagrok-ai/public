import {defineConfig, devices} from '@playwright/test';

export default defineConfig({
  testDir: '.',
  testMatch: '**/*-spec.ts',
  // helpers/session-spec.ts requires manual password entry in a visible
  // browser; it would hang any headless run. Opt-in via HEADED_MANUAL=1 +
  // explicit file path (the in-file `test.skip` handles the runtime guard).
  testIgnore: ['**/helpers/session-spec.ts'],
  // Run all eligible tests in parallel — within and across files. Specs that
  // share UI state across multiple test() blocks should opt out per-file
  // with `test.describe.configure({mode: 'serial'})`.
  fullyParallel: true,
  // 4 CI workers (down from 8): each worker drives its own Chromium hitting
  // the same single-replica CI Datlas, and at 8× concurrency the
  // `.grok-preloader` wait in spec-login.ts:33 was timing out under cold-start
  // contention (34/154 specs). 4 is the empirical sweet spot.
  workers: process.env.CI ? 4 : 1,
  retries: process.env.CI ? 1 : 0,
  // Per-test default timeout (1 minute). Specs that need longer override
  // explicitly via `test.setTimeout(N)` at the top of the test body.
  timeout: 60_000,
  expect: {timeout: 15_000},
  reporter: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME
    ? [['list'], ['json', {outputFile: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME}]]
    : [['list']],
  outputDir: '../../test-playwright-output',
  use: {
    baseURL: process.env.DATAGROK_URL ?? 'http://localhost:8888',
    viewport: {width: 1920, height: 1080},
    actionTimeout: 15_000,
    navigationTimeout: 60_000,
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
  },
  projects: [
    {name: 'chromium', use: {...devices['Desktop Chrome']}},
  ],
});
