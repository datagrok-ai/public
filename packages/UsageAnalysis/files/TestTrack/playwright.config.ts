import {defineConfig, devices} from '@playwright/test';

export default defineConfig({
  testDir: '.',
  testMatch: '**/*-spec.ts',
  // fullyParallel stays off so tests *within* a file remain sequential (most
  // specs share UI state — e.g. a project they create in step 1 and assert on
  // in step 2). Spec *files* still run in parallel across `workers`.
  fullyParallel: false,
  workers: process.env.CI ? 4 : 1,
  retries: process.env.CI ? 1 : 0,
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
