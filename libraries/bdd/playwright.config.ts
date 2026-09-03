/* The one Playwright config every bdd project runs with (`grok-bdd run` passes it as --config and
   points BDD_ROOT at the project). The shared Datagrok base (login storage state, viewport,
   traces) with the project's generated/ as the test dir. */
import {existsSync} from 'node:fs';
import {basename, dirname, join, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';
import {defineConfig} from '@playwright/test';
import {baseConfig} from '@datagrok-libraries/test/src/playwright/base-config.js';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(process.env.BDD_ROOT ?? (basename(here) === 'dist' ? dirname(here) : here));
const url = (process.env.DATAGROK_URL ?? 'http://localhost:8888').replace(/\/$/, '');
const globalSetup = ['global-setup.js', 'global-setup.ts'].map((f) => join(here, 'src', 'runtime', f)).find(existsSync)!;

export default defineConfig({
  ...baseConfig,
  testDir: join(root, 'generated'),
  outputDir: join(root, 'test-results'),
  globalSetup,
  use: {
    ...baseConfig.use,
    baseURL: url,
    storageState: join(root, 'e2e', '.auth.json'),
    launchOptions: {args: [`--unsafely-treat-insecure-origin-as-secure=${url}`]},
  },
  projects: [{name: 'bdd'}],
});
