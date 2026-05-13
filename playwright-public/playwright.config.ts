import {defineConfig, devices} from '@playwright/test';

export default defineConfig({
  testDir: '.',
  testMatch: '**/*.test.ts',
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
  },
  projects: [
    {name: 'chromium', use: {...devices['Desktop Chrome']}},
  ],
});
