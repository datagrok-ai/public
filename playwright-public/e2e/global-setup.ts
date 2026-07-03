import {chromium, FullConfig} from '@playwright/test';
import * as fs from 'fs';
import * as path from 'path';

// Global setup: turns the `DATAGROK_AUTH_TOKEN` provided by `grok test` into
// browser storage state (cookie + localStorage) and writes it to every path
// the existing specs reference, so each suite picks it up without changes.
//
// `DATAGROK_AUTH_TOKEN` and `DATAGROK_URL` are set by `runPlaywrightTests` in
// `public/tools/bin/utils/playwright-runner.ts` before spawning Playwright;
// we throw early if either is missing.
export default async function globalSetup(_config: FullConfig) {
  const baseURL = process.env.DATAGROK_URL;
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!baseURL || baseURL.length === 0)
    throw new Error('DATAGROK_URL is not set. Run via `grok test --skip-puppeteer` from this directory.');
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test --skip-puppeteer` — it derives the token from ~/.grok/config.yaml.');

  const browser = await chromium.launch();
  try {
    const ctx = await browser.newContext();
    const page = await ctx.newPage();

    // Mirror what `grok test`'s Puppeteer pass does (test-utils.ts:135):
    // navigate to /oauth/ first so cookies attach to the right origin, then
    // drop the auth token into both cookie and localStorage.
    //
    // `waitUntil: 'domcontentloaded'` + a generous timeout: the CI client is a
    // cold `pub serve` (debug Dart) behind nginx that compiles bundles on first
    // request, so the `load` event (all dart.js fetched) can exceed Playwright's
    // default 30 s — this is what timed out global-setup on the CI client while
    // it succeeds on warm dev. We only need the document parsed here; the real
    // readiness gate is the preloader + Browse waits after the reload below.
    const gotoOpts = {waitUntil: 'domcontentloaded' as const, timeout: 120_000};
    await page.goto(baseURL.replace(/\/$/, '') + '/oauth/', gotoOpts);
    const u = new URL(baseURL);
    await ctx.addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
    await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
    await page.goto(baseURL, gotoOpts);
    await page.waitForFunction(() => document.querySelector('#grok-preloader, .grok-preloader') == null, undefined, {timeout: 30_000}).catch(() => {});
    await page.addStyleTag({content: '#grok-preloader,.grok-preloader{pointer-events:none!important}.d4-tooltip{display:none!important}'}).catch(() => {});
    await page.locator('[name="Browse"]').first().waitFor({timeout: 60_000});

    const state = await ctx.storageState();
    const root = path.resolve(__dirname, '..');
    const e2eDir = path.join(root, 'e2e');
    if (!fs.existsSync(e2eDir))
      fs.mkdirSync(e2eDir, {recursive: true});

    // Three different paths are referenced by the existing specs:
    //   connections/* → AUTH_STATE = 'e2e/.auth.json'        (relative to CWD)
    //   queries/*     → AUTH_STATE = 'e2e/.auth.public.json' (relative to CWD)
    //   scripts/*     → path.resolve(__dirname, '..', '.auth.json')
    // CWD when Playwright runs under `grok test` is this package directory,
    // so all three resolve to siblings of this file.
    const json = JSON.stringify(state);
    fs.writeFileSync(path.join(e2eDir, '.auth.json'), json);
    fs.writeFileSync(path.join(e2eDir, '.auth.public.json'), json);
    fs.writeFileSync(path.join(root, '.auth.json'), json);
  }
  finally {
    await browser.close();
  }
}
