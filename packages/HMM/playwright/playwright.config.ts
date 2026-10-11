import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// HMM package-owned E2E suite. General config lives in the shared base; only run-dir
// specifics belong here. The screenshot tour (tour/*.tour.ts, explained screenshots of
// every feature) runs only on request:
//   HMM_TOUR=1 grok test --host localhost --skip-build --skip-puppeteer --category tour
export default defineConfig({
  ...baseConfig,
  testDir: '.',
  testMatch: process.env.HMM_TOUR ? '**/*.tour.ts' : '**/*.test.ts',
});
