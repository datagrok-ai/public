import {defineConfig, devices} from '@playwright/test';

export default defineConfig({
  testDir: '.',
  testMatch: '**/*.test.ts',
  // `helpers/` hosts utility modules + the headed-only `session-helpers.test.ts`
  // that validates `logoutAndLoginAs` against a real second user. It cannot
  // run headlessly (manual password entry) and is not a project-feature test,
  // so we exclude the whole folder from test discovery.
  testIgnore: ['**/helpers/**'],
  // Many specs share UI/server state across tests in a file (connection lifecycle,
  // query lifecycle, scripts CRUD). Keep one worker by default — the suite is
  // designed to be sequential. CI can opt into parallelism per-file via
  // `test.describe.parallel` if/when specs are made independent.
  fullyParallel: false,
  workers: 1,
  retries: process.env.CI ? 1 : 0,
  // Per-test default timeout (2 minutes). Some Browse-tree drill-downs and
  // identifier-config flows wait on cold-cache server fetches.
  timeout: 120_000,
  expect: {timeout: 15_000},
  globalSetup: './e2e/global-setup.ts',
  reporter: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME
    ? [['list'], ['json', {outputFile: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME}]]
    : [['list']],
  outputDir: 'test-output',
  use: {
    baseURL: process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai',
    storageState: 'e2e/.auth.json',
    viewport: {width: 1920, height: 1080},
    actionTimeout: 15_000,
    navigationTimeout: 60_000,
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
    // Grant clipboard access so headless-CI copy flows (copy-as-HELM, sketcher
    // Copy as SMILES/MOLBLOCK) can use navigator.clipboard instead of throwing.
    permissions: ['clipboard-read', 'clipboard-write'],
    launchOptions: {
      // The CI stack serves over plain HTTP (xamgle-nginx:8889). navigator.clipboard
      // only exists in a secure context, so copy flows fail with "writeText on null".
      // Treat the target origin as secure to expose the Clipboard API on HTTP.
      args: [`--unsafely-treat-insecure-origin-as-secure=${(process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai').replace(/\/$/, '')}`],
    },
  },
  projects: [
    {name: 'chromium', use: {...devices['Desktop Chrome']}},
  ],
});
