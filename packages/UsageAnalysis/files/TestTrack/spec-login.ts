import {test, Page} from '@playwright/test';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import * as yaml from 'js-yaml';

export const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';

export const specTestOptions = {
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
};

export interface StepError { step: string; error: string; }

export const stepErrors: StepError[] = [];

export async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

// Wait until the top-menu Chem entry registers — it appears only after the
// Molecule semType is detected on the active table and the Chem package is
// ready. Polls up to ~15s for the `[name="div-Chem"]` element to be present.
export async function waitForChemMenu(page: Page) {
  await page.locator('[name="div-Chem"]').first().waitFor({state: 'attached', timeout: 15_000});
}

// Poll until a column with Molecule semType exists on the active table (or the
// spec's window.__df handle). The Chem autostart detector runs asynchronously
// AFTER the Chem menu attaches, so waitForChemMenu alone does not guarantee
// semType has been applied — checking immediately races the detector.
export async function waitForMolecule(page: Page, timeoutMs = 45_000) {
  await page.waitForFunction(() => {
    const g = (window as any).grok;
    const tables = [g?.shell?.t, (window as any).__df].filter(Boolean);
    return tables.some((t: any) => t.columns.toList().some((c: any) => c.semType === 'Molecule'));
  }, null, {timeout: timeoutMs});
}

async function injectToken(page: Page, token: string) {
  // Navigate to the origin first so the cookie/localStorage entries are
  // attached to the right host. The `/oauth/` path matches what `grok test`
  // does for Puppeteer (test-utils.ts:135).
  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl);
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

export async function loginToDatagrok(page: Page) {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  await injectToken(page, token);
}

// Build the direct file-browse URL for the configured instance (DATAGROK_URL —
// dev/public/release/local). `relPath` is dot-namespaced, e.g.
// 'System.AppData/Helm/samples/helm-showcase.csv'.
export function fileBrowseUrl(relPath: string): string {
  return `${baseUrl}/file/${relPath}?browse=files`;
}

// Authenticate, then navigate DIRECTLY to the dataset's file-browse URL so the
// platform and the dataset open in a SINGLE navigation — no open-platform-then-
// readCsv round-trip. The URL is derived from the configured instance, so the
// same spec opens the right dataset on dev / public / release without changes.
// Waits for the preloader to clear and the grid viewer to attach.
export async function loginAndOpenFile(page: Page, relPath: string) {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  // Set the auth cookie/localStorage at the origin first…
  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  // …then go straight to the dataset URL (platform + dataset in one navigation).
  await page.goto(fileBrowseUrl(relPath));
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60_000});
}

// Read a second-user dev key (`key2:`) from ~/.grok/config.yaml for the server
// whose url matches the current DATAGROK_URL (falling back to the configured
// default server). This makes the second-user login work under a plain
// `grok test` even when the installed runner doesn't forward key2 — the only
// thing the installed runner lacks is this config fallback.
function readSecondUserDevKeyFromConfig(): {apiUrl: string; key2: string} | null {
  try {
    const confPath = path.join(os.homedir(), '.grok', 'config.yaml');
    if (!fs.existsSync(confPath)) return null;
    const cfg = yaml.load(fs.readFileSync(confPath, 'utf8')) as any;
    const servers = cfg?.servers ?? {};
    let wantHost: string | null = null;
    try { wantHost = new URL(baseUrl).host; } catch (_) { wantHost = null; }
    for (const name of Object.keys(servers)) {
      const s = servers[name];
      if (!s?.url || !s?.key2) continue;
      let h: string | null = null;
      try { h = new URL(s.url).host; } catch (_) { h = null; }
      if (h && wantHost && h === wantHost)
        return {apiUrl: String(s.url).replace(/\/$/, ''), key2: String(s.key2)};
    }
    const def = cfg?.default;
    if (def && servers[def]?.key2 && servers[def]?.url)
      return {apiUrl: String(servers[def].url).replace(/\/$/, ''), key2: String(servers[def].key2)};
    return null;
  } catch (_) {
    return null;
  }
}

async function exchangeDevKeyForToken(apiUrl: string, key: string): Promise<string> {
  const resp = await fetch(`${apiUrl}/users/login/dev/${key}`, {method: 'POST'});
  const json = await resp.json() as any;
  if (json?.isSuccess === true && json?.token) return json.token;
  throw new Error(`Second-user dev-key login failed at ${apiUrl}: ${JSON.stringify(json).slice(0, 200)}`);
}

// Resolve the second-user token: env first (set by the runner), else exchange
// the config `key2:` dev key for a token. Throws when neither is available —
// a two-user spec MUST NOT silently pass without its second user. Cached so
// the login claim can be read (getSecondUserLogin) without a second exchange.
let _secondTokenCache: string | null = null;
export async function resolveSecondUserToken(): Promise<string> {
  if (_secondTokenCache) return _secondTokenCache;
  const envTok = process.env.DATAGROK_AUTH_TOKEN_2;
  if (envTok && envTok.length > 0) return (_secondTokenCache = envTok);
  const cfg = readSecondUserDevKeyFromConfig();
  if (!cfg)
    throw new Error(
      'No second-user credentials available. Set DATAGROK_AUTH_TOKEN_2 / DATAGROK_DEV_KEY_2, ' +
      'or add a `key2:` (second-user dev key) to the matching server in ~/.grok/config.yaml.');
  return (_secondTokenCache = await exchangeDevKeyForToken(cfg.apiUrl, cfg.key2));
}

// Read the second user's login from the JWT `sub` (/`usr.login`) claim — no
// page switch needed. Lets a two-user spec learn WHO the recipient is so it
// can share with them, without the costly "switch in, read, switch back"
// probe (which added two page reloads per spec).
export async function getSecondUserLogin(): Promise<string> {
  const token = await resolveSecondUserToken();
  try {
    const payload = token.replace(/^Bearer\s+/i, '').split('.')[1];
    const claims = JSON.parse(Buffer.from(payload, 'base64').toString('utf8'));
    const login = claims?.sub ?? claims?.usr?.login;
    if (!login) throw new Error('no sub/usr.login claim');
    return login;
  } catch (e: any) {
    throw new Error(`Could not read second-user login from token claim: ${e?.message ?? e}`);
  }
}

export async function loginAsSecondUser(page: Page) {
  const token2 = await resolveSecondUserToken();
  await injectToken(page, token2);
}
