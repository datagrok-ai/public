import {test, Page} from '@playwright/test';

export const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';

export const specTestOptions = {
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
};

export interface StepError { step: string; error: string; }

export const stepErrors: StepError[] = [];

export async function softStep(name: string, fn: () => Promise<void>): Promise<void> {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

// Wait until the top-menu Chem entry registers — it appears only after the
// Molecule semType is detected on the active table and the Chem package is
// ready. Polls up to ~15s for the `[name="div-Chem"]` element to be present.
export async function waitForChemMenu(page: Page): Promise<void> {
  await page.locator('[name="div-Chem"]').first().waitFor({state: 'attached', timeout: 15_000});
}

// Poll until a column with Molecule semType exists on the active table (or the
// spec's window.__df handle). The Chem autostart detector runs asynchronously
// AFTER the Chem menu attaches, so waitForChemMenu alone does not guarantee
// semType has been applied — checking immediately races the detector.
export async function waitForMolecule(page: Page, timeoutMs = 45_000): Promise<void> {
  await page.waitForFunction(() => {
    const g = (window as any).grok;
    const tables = [g?.shell?.t, (window as any).__df].filter(Boolean);
    return tables.some((t: any) => t.columns.toList().some((c: any) => c.semType === 'Molecule'));
  }, null, {timeout: timeoutMs});
}

async function injectToken(page: Page, token: string): Promise<void> {
  // Navigate to the origin first so the cookie/localStorage entries are
  // attached to the right host. The `/oauth/` path matches what `grok test`
  // does for Puppeteer (test-utils.ts:135).
  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl);
  // Cold CI Datlas keeps #grok-preloader in the DOM past the timeout even though
  // the Browse tree is already interactive (confirmed via CI page snapshot). Make
  // the wait best-effort and neutralise the preloader's click-interception so a
  // lingering preloader can neither hard-fail login nor intercept later clicks.
  await page.waitForFunction(() => document.querySelector('#grok-preloader, .grok-preloader') == null, null, {timeout: 30_000})
    .catch(() => { /* tolerate a lingering preloader — neutralised below */ });
  await page.addStyleTag({content: `
    #grok-preloader, .grok-preloader { pointer-events: none !important; }
    .d4-tooltip { display: none !important; }
  `}).catch(() => {});
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

export async function loginToDatagrok(page: Page): Promise<void> {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  await injectToken(page, token);
}

// Direct file-browse URL for a dataset (dot-namespaced relPath, e.g.
// System.AppData/Helm/samples/HELM.csv). Ported 1:1 from TestTrack spec-login.
export function fileBrowseUrl(relPath: string): string {
  return `${baseUrl}/file/${relPath}?browse=files`;
}

// Authenticate, then open the dataset into a table view. The TestTrack original
// deep-linked to the dot-namespaced file-browse URL (/file/System.AppData/...?browse=files),
// but that route 404s on the minimal ui_tests CI stack (the xamgle SPA resolves it as a
// static asset). Instead load the platform normally and open the file via the file API —
// the CI-safe pattern used across the suite (cf. Bio openBioDataset). The dot-namespaced
// connector prefix (System.AppData) is converted to the API colon form (System:AppData).
export async function loginAndOpenFile(page: Page, relPath: string): Promise<void> {
  await loginToDatagrok(page);
  const apiPath = relPath.replace(/^System\.AppData\//, 'System:AppData/');
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    const df = await g.dapi.files.readCsv(p);
    g.shell.addTableView(df);
  }, apiPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60_000});
}

// Resolve the second-user token: env first (DATAGROK_AUTH_TOKEN_2, which the CI runner exports after
// provisioning the `test2` user), else exchange a second-user dev key (DATAGROK_DEV_KEY_2) for a token.
// Throws when neither is available — a two-user spec MUST NOT silently pass without its second user.
// Cached so the login claim can be read (getSecondUserLogin) without a second exchange.
let _secondTokenCache: string | null = null;
export async function resolveSecondUserToken(): Promise<string> {
  if (_secondTokenCache) return _secondTokenCache;
  const envTok = process.env.DATAGROK_AUTH_TOKEN_2;
  if (envTok && envTok.length > 0) return (_secondTokenCache = envTok);
  const key2 = process.env.DATAGROK_DEV_KEY_2;
  if (key2 && key2.length > 0) {
    const apiUrl = (process.env.DATAGROK_URL ?? baseUrl).replace(/\/$/, '') + '/api';
    const resp = await fetch(`${apiUrl}/users/login/dev/${key2}`, {method: 'POST'});
    const json = await resp.json() as any;
    if (json?.isSuccess === true && json?.token) return (_secondTokenCache = json.token);
    throw new Error(`Second-user dev-key login failed at ${apiUrl}: ${JSON.stringify(json).slice(0, 200)}`);
  }
  throw new Error('No second-user credentials available. Set DATAGROK_AUTH_TOKEN_2 (the CI runner provides it) or DATAGROK_DEV_KEY_2.');
}

// Resolve the second user's login so a two-user spec can learn WHO the recipient is in order to share
// the project with them. Prefer decoding the JWT `sub` (/`usr.login`) claim (no network round-trip);
// fall back to /api/users/current with the token when it isn't a decodable JWT (e.g. a plain session
// token from the CI runner's /users/login provisioning of `test2`).
export async function getSecondUserLogin(): Promise<string> {
  const token = await resolveSecondUserToken();
  try {
    const payload = token.replace(/^Bearer\s+/i, '').split('.')[1];
    if (payload) {
      const claims = JSON.parse(Buffer.from(payload, 'base64').toString('utf8'));
      const login = claims?.sub ?? claims?.usr?.login;
      if (login) return login;
    }
  } catch (_) { /* not a JWT — fall through to the REST lookup */ }
  const apiUrl = (process.env.DATAGROK_URL ?? baseUrl).replace(/\/$/, '') + '/api';
  const resp = await fetch(`${apiUrl}/users/current`, {headers: {Authorization: token}});
  const user = await resp.json() as any;
  if (!user?.login)
    throw new Error('Could not resolve second-user login (JWT claim and /users/current both failed)');
  return user.login;
}

export async function loginAsSecondUser(page: Page): Promise<void> {
  const token2 = await resolveSecondUserToken();
  await injectToken(page, token2);
}
