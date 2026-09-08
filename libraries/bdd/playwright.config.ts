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
    // a failed run keeps its trace (actions, console, network) and the failure screenshot, but
    // neither DOM snapshots (serializing the shell's DOM around every action was ~45% of a
    // feature's time) nor a screenshot per action (~12 s over six features, only a filmstrip);
    // `grok-bdd run --trace on` records everything, `--video on` a video
    trace: {mode: 'retain-on-failure', snapshots: false, screenshots: false},
    // a GPU-rasterized canvas re-rasterizes on the CPU after enough pixel readbacks, and that first
    // paint differs in antialiasing from the one before it (~2000 px on a bar chart) — headed only,
    // headless is software-rasterized throughout
    launchOptions: {args: [`--unsafely-treat-insecure-origin-as-secure=${url}`, '--disable-accelerated-2d-canvas']},
  },
  projects: [{name: 'bdd'}],
});
