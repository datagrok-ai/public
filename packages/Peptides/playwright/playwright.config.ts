import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config';

// Peptides package-owned E2E suite. General config lives in the shared base; only run-dir
// specifics belong here.
export default defineConfig({
  ...baseConfig,
  testDir: '.',
});
