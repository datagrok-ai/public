import {defineConfig, devices} from '@playwright/test';
import * as path from 'path';

// Shared Playwright config for every Datagrok E2E suite — the playwright-public core
// suites and the package-owned `playwright/` folders alike. A consumer config does:
//
//   import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';
//   export default defineConfig({...baseConfig, testDir: '.'});
//
// Put only genuinely package-specific overrides in the consumer config; everything
// general lives here.
const DATAGROK_URL = (process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai').replace(/\/$/, '');

export const baseConfig = defineConfig({
  testMatch: '**/*.test.ts',
  // Many specs share UI/server state across tests in a file, so the suite is
  // designed to be sequential. Consumers can opt into parallelism per-file via
  // test.describe.parallel once their specs are independent.
  fullyParallel: false,
  workers: 1,
  retries: process.env.CI ? 1 : 0,
  timeout: 120_000,
  expect: {timeout: 15_000},
  // Absolute path so it resolves to this lib regardless of which consumer config
  // imports baseConfig (Playwright would otherwise resolve a relative globalSetup
  // against the consumer's config dir).
  globalSetup: path.join(__dirname, 'global-setup'),
  reporter: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME
    ? [['list'], ['json', {outputFile: process.env.PLAYWRIGHT_JSON_OUTPUT_NAME}]]
    : [['list']],
  outputDir: 'test-output',
  use: {
    baseURL: DATAGROK_URL,
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
      // The CI stack serves over plain HTTP; navigator.clipboard only exists in a
      // secure context. Treat the target origin as secure to expose the Clipboard API.
      args: [`--unsafely-treat-insecure-origin-as-secure=${DATAGROK_URL}`],
    },
  },
  projects: [
    {name: 'chromium', use: {...devices['Desktop Chrome']}},
  ],
});
