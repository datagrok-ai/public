import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// Charts package-owned E2E suite. All general config lives in the shared base
// (@datagrok-libraries/test/src/playwright/base-config); only run-dir specifics here.
export default defineConfig({
  ...baseConfig,
  testDir: '.',
});
