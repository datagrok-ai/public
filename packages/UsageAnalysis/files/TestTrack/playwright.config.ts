import {defineConfig, devices} from '@playwright/test';

export default defineConfig({
  testDir: '.',
  testMatch: '**/*-spec.ts',
  // Run all eligible tests in parallel — within and across files. Specs that
  // share UI state across multiple test() blocks should opt out per-file
  // with `test.describe.configure({mode: 'serial'})`.
  fullyParallel: true,
  workers: process.env.CI ? 8 : 1,
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
