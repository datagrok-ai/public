import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// HMM package-owned E2E suite. General config lives in the shared base; only run-dir
// specifics belong here. The screenshot tour (tour/*.tour.ts, explained screenshots of
// every feature) runs only on request:
//   HMM_TOUR=1 grok test --host localhost --skip-build --skip-puppeteer --category tour
// and so does the performance run at scale (perf/*.perf.ts, 100,000 generated antibodies):
//   HMM_PERF=1 grok test --host localhost --skip-build --skip-puppeteer --category perf
export default defineConfig({
  ...baseConfig,
  testDir: '.',
  testMatch: process.env.HMM_PERF ? '**/*.perf.ts' : process.env.HMM_TOUR ? '**/*.tour.ts' : '**/*.test.ts',
});
