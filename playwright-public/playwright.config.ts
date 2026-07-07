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
});
