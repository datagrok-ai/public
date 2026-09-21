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
// Falling back to a deployed environment is not a safe default: specs create projects,
// edit shares and delete entities, so a harness that forgot to export DATAGROK_URL would
// do all of that against dev while the output still looked like a normal run. `grok test`
// always sets it; anything else should say where it is pointing.
const DATAGROK_URL = (process.env.DATAGROK_URL ?? 'http://localhost:8888').replace(/\/$/, '');
// The CI stack is plain HTTP, so an insecure origin has to be allow-listed as secure or the
// whole secure-context surface is missing. Measured on Chrome-for-Testing 153 (playwright
// chromium v1243): the default `chromium-headless-shell` IGNORES
// --unsafely-treat-insecure-origin-as-secure (isSecureContext stays false), while the full
// Chrome build honours it headless. Every spec takes that build: handing the shell to the
// specs that looked like they only needed painting cost the Chem substructure query and the
// sketcher paste, which both sit on that surface.
const INSECURE_TARGET = DATAGROK_URL.startsWith('http://');

// `test.use({launchOptions})` REPLACES this block rather than merging into it, so every
// consumer that sets launch options of its own (specTestOptions) must spread this in.
export const secureOriginLaunchOptions = {
  ...(INSECURE_TARGET ? {channel: 'chromium'} : {}),
  args: [`--unsafely-treat-insecure-origin-as-secure=${DATAGROK_URL}`],
};

/** Kept as the name the clipboard specs opt into; the build is the same one everything takes. */
export const clipboardLaunchOptions = secureOriginLaunchOptions;

export const baseConfig = defineConfig({
  testMatch: '**/*.test.ts',
  // Many specs share UI/server state across tests in a file, so ordering WITHIN a file
  // stays sequential. Separate files are independent, so they can run alongside each
  // other; PLAYWRIGHT_WORKERS lets CI dial that per agent without a republish.
  fullyParallel: false,
  workers: Number(process.env.PLAYWRIGHT_WORKERS ?? 4),
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
    launchOptions: secureOriginLaunchOptions,
  },
  projects: [
    {name: 'chromium', use: {...devices['Desktop Chrome']}},
  ],
});
